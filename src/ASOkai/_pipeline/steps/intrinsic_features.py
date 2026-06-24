#!/usr/bin/env python
"""
Filename: src/ASOkai/_pipeline/steps/intrinsic_features.py
Author: Arash Ayat
Copyright: 2026, Alexander Schliep
Version: 0.1.1
Description: Definition and CLI entrypoint for the intrinsic-features step.
License: LGPL-3.0-or-later
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from importlib.resources import files

from ASOkai.Analysis import IntrinsicFeaturesAnalysis
from ASOkai._pipeline.base import AnalysisStep


class IntrinsicFeaturesStep(AnalysisStep):
    name = "intrinsic-features"
    description = "Computes intrinsic features (GC content, T-runs, AT-runs) for each ASO target site."
    analysis_cls = IntrinsicFeaturesAnalysis
    cli_module = "ASOkai._pipeline.steps.intrinsic_features"
    dependencies: list[str] = ["create-target-gene"]
    config_map = {
        "target_id":   "target.target_id",
        "target_name": "target.target_name",
        "k":           "target.k",
        "region":      "target.region",
        "assembly":    "genome.assembly_id",
    }
    input_overrides: dict[str, str] = {
        "target_gene": "target.target_gene_path",
    }

    @property
    def cwl_path(self) -> str:
        return str(files("ASOkai._cwl.steps").joinpath("intrinsic-features.cwl"))

    def _effective_target_id(self, config: dict) -> str:
        return config["target"].get("target_id") or config["target"].get("target_name")

    def _analysis_dir(self, config: dict) -> Path:
        assembly = config["genome"]["assembly_id"]
        target_id = self._effective_target_id(config)
        return (
            Path(config["datadir"])
            / assembly
            / "targets"
            / target_id
            / "analysis"
            / "intrinsic"
        )

    def outdir(self, config: dict) -> Path:
        return self._analysis_dir(config)

    def output_paths(self, config: dict) -> dict[str, Path]:
        base = self._analysis_dir(config)
        target_id = self._effective_target_id(config)
        k = config["target"]["k"]
        region = config["target"]["region"]
        return {
            "intrinsic_features": base / f"{target_id}_k{k}_{region}_intrinsic.json",
        }

    def outputs_exist(self, config: dict) -> bool:
        return all(p.exists() for p in self.output_paths(config).values())

    def cleanup(self, config: dict) -> None:
        for p in self.output_paths(config).values():
            if p.exists():
                p.unlink()

    def load_analysis_inputs(self, args) -> dict:
        from ASOkai.Targets import TargetGene

        return {
            "target_gene": TargetGene.from_file(str(args.target_gene)),
        }

    def analysis_kwargs(self, args, inputs: dict) -> dict:
        return {
            "sites": inputs["target_gene"].sites,
        }

    def analysis_metadata(self, args, inputs: dict) -> dict:
        return {
            "analysis": self.name,
            "assembly": args.assembly,
            "target_id": args.target_id,
            "target_name": args.target_name,
            "k": args.k,
            "region": args.region,
        }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute intrinsic features for ASO target sites.",
    )
    parser.add_argument("--target-gene", required=True, type=Path, help="Path to target gene JSON.")
    parser.add_argument("--assembly", required=True, help="Assembly ID (e.g. GRCh38).")
    parser.add_argument("--target-id", default=None, help="Ensembl gene ID. Takes priority over --target-name.")
    parser.add_argument("--target-name", default=None, help="Gene name. Used if --target-id is not provided.")
    parser.add_argument("--k", required=True, type=int, help="ASO length.")
    parser.add_argument("--region", required=True, choices=["exonic_only", "pre-mrna", "transcriptomic"], help="Target region type.")
    parser.add_argument("--output", required=True, type=Path, help="Full path for the output JSON file.")
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint called by CWL baseCommand: intrinsic-features."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    if not args.target_id and not args.target_name:
        parser.error("Either --target-id or --target-name is required.")

    return IntrinsicFeaturesStep().run_from_args(args)


if __name__ == "__main__":
    sys.exit(main())
