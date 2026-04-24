"""Tests for the public ASOkai CLI entrypoint."""

from click.testing import CliRunner

from ASOkai.cli.main import main


def test_hidden_step_command_is_not_listed_in_help():
    runner = CliRunner()

    result = runner.invoke(main, ["--help"])

    assert result.exit_code == 0
    assert " step " not in result.output


def test_step_command_dispatches_to_builtin_step(monkeypatch):
    from pipeline.steps import download_genome

    received: dict[str, list[str] | None] = {}

    def fake_main(argv: list[str] | None = None) -> int:
        received["argv"] = argv
        return 0

    monkeypatch.setattr(download_genome, "main", fake_main)

    runner = CliRunner()
    result = runner.invoke(
        main,
        ["step", "download-genome", "--assembly", "GRCh38", "--release", "114"],
    )

    assert result.exit_code == 0
    assert received["argv"] == ["--assembly", "GRCh38", "--release", "114"]


def test_step_command_rejects_unknown_steps():
    runner = CliRunner()

    result = runner.invoke(main, ["step", "not-a-step"])

    assert result.exit_code != 0
    assert "Unknown step 'not-a-step'" in result.output
