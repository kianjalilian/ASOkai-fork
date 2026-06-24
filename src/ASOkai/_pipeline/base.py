#!/usr/bin/env python
"""
Filename: src/ASOkai/_pipeline/base.py
Author: Arash Ayat
Copyright: 2026, Alexander Schliep
Version: 0.1.1
Description: Base protocols that steps, tasks, and workflows must implement.
License: LGPL-3.0-or-later
"""
from __future__ import annotations

from abc import ABC, abstractmethod
import json
from pathlib import Path
from typing import Any, ClassVar, Protocol, runtime_checkable


@runtime_checkable
class Runnable(Protocol):
    """
    Shared contract for every named pipeline unit (step, task, or workflow).
    """

    name: str
    description: str

    def output_paths(self, config: dict) -> dict[str, Path]:
        ...

    def outputs_exist(self, config: dict) -> bool:
        ...

    def cleanup(self, config: dict) -> None:
        ...


class Step(Runnable, ABC):
    """Atomic pipeline unit backed by a CWL command-line tool."""

    name: ClassVar[str]
    description: ClassVar[str]
    dependencies: ClassVar[list[str]]
    config_map: ClassVar[dict[str, str]]
    input_overrides: ClassVar[dict[str, str]]
    cli_module: ClassVar[str]

    @property
    @abstractmethod
    def cwl_path(self) -> str:
        """Path to the static CWL file for this step."""
        ...

    @abstractmethod
    def outdir(self, config: dict) -> Path:
        """Directory passed to the CWL executor for this step's outputs."""
        ...

    @abstractmethod
    def output_paths(self, config: dict) -> dict[str, Path]:
        ...

    @abstractmethod
    def outputs_exist(self, config: dict) -> bool:
        ...

    @abstractmethod
    def cleanup(self, config: dict) -> None:
        ...


class CoreStep(Step):
    """Pipeline step that prepares, downloads, builds, or transforms core data."""


class AnalysisStep(Step):
    """Pipeline step that runs analysis logic and writes analysis results."""

    analysis_cls: ClassVar[type | None] = None

    def load_analysis_inputs(self, args) -> dict[str, Any]:
        """Load input objects needed to construct the analysis."""
        return {}

    def analysis_kwargs(self, args, inputs: dict[str, Any]) -> dict[str, Any]:
        """Build keyword arguments for the configured analysis class."""
        return {}

    def analysis_metadata(self, args, inputs: dict[str, Any]) -> dict[str, Any]:
        """Build output metadata written next to analysis results."""
        return {"analysis": self.name}

    def write_analysis_output(self, args, payload: dict[str, Any]) -> None:
        """Write the analysis payload."""
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=4))

    def run_from_args(self, args) -> int:
        """Load inputs, run the configured analysis, and write its JSON output."""
        if self.analysis_cls is None:
            raise RuntimeError(f"Analysis step '{self.name}' does not define analysis_cls.")

        inputs = self.load_analysis_inputs(args)
        analysis = self.analysis_cls(**self.analysis_kwargs(args, inputs))
        payload = {
            **self.analysis_metadata(args, inputs),
            "results": analysis.run(),
        }
        self.write_analysis_output(args, payload)
        return 0


@runtime_checkable
class Task(Runnable, Protocol):
    """Ordered composition of Steps."""

    steps: list[Step]


@runtime_checkable
class Workflow(Runnable, Protocol):
    """
    Ordered composition of Runnables (Steps, Tasks, or nested Workflows).

    CWL is generated at runtime by recursively flattening ``members`` to Steps.
    """

    members: list[Runnable]
