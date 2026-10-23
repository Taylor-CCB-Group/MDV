from importlib import metadata

import click
import pytest
from click.testing import CliRunner

import mdvtools.cli as cli_module
from mdvtools.diagnostics import (
    format_error_report,
    format_version_info,
    get_version_info,
    report_errors,
)


def test_version_flag_prints_version_and_exits_cleanly():
    runner = CliRunner()
    for flag in ("--version", "-V"):
        result = runner.invoke(cli_module.cli, [flag])
        assert result.exit_code == 0, result.output
        assert result.output.startswith(f"mdvtools {metadata.version('mdvtools')}")
        assert "python" in result.output
        assert "platform" in result.output


def test_format_version_info_includes_build_and_dependencies():
    info = get_version_info(include_dependencies=True)
    info["build"] = {
        "git_commit_hash": "0123456789abcdef",
        "git_branch": "main",
        "git_commit_date": "2026-01-01",
        "git_dirty": True,
    }
    text = format_version_info(info)
    assert "commit    0123456789ab (main, 2026-01-01, uncommitted changes)" in text
    assert f"click {metadata.version('click')}" in text


def test_unexpected_error_prints_report_and_exits_1(monkeypatch):
    runner = CliRunner()
    result = runner.invoke(cli_module.cli, ["serve", "/definitely/not/a/project"])
    assert result.exit_code == 1
    # the report replaces click's default handling, so no raw exception escapes
    assert isinstance(result.exception, SystemExit)
    out = result.output
    assert "mdvtools error report" in out
    assert "Traceback (most recent call last)" in out
    assert "in serve" in out
    assert "in wrapper" not in out
    assert f"mdvtools {metadata.version('mdvtools')}" in out
    assert out.rstrip().splitlines()[-2] == (
        "Error: FileNotFoundError: /definitely/not/a/project not found"
    )


def test_usage_errors_are_left_to_click():
    runner = CliRunner()
    result = runner.invoke(cli_module.cli, ["serve"])
    assert result.exit_code == 2
    assert "mdvtools error report" not in result.output
    assert "Missing argument" in result.output


def test_report_errors_passes_through_click_exceptions():
    @report_errors
    def fails():
        raise click.BadParameter("nope")

    with pytest.raises(click.BadParameter):
        fails()


def test_format_error_report_survives_version_info_failure(monkeypatch):
    import mdvtools.diagnostics as diagnostics

    def broken(**_kwargs):
        raise RuntimeError("no metadata")

    monkeypatch.setattr(diagnostics, "get_version_info", broken)
    try:
        raise ValueError("boom")
    except ValueError as exc:
        report = format_error_report(exc, argv=["/usr/bin/mdvtools", "serve", "my project"])
    assert "command:  mdvtools serve 'my project'" in report
    assert "version information unavailable" in report
    assert "Error: ValueError: boom" in report
