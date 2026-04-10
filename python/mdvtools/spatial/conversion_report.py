from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import traceback
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path


POINT_TRANSFORM_CHOICES = [
    "image",
    "auto",
    "xenium",
    "identity",
    "annotated-element",
]


@dataclass
class ConversionResult:
    dataset_name: str
    dataset_path: str
    status: str
    returncode: int
    duration_seconds: float
    command: list[str]
    stdout_log: str
    stderr_log: str
    warning_count: int
    warning_samples: list[str]
    error_type: str | None
    error_message: str | None
    error_signature: str | None
    traceback_excerpt: str | None


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run mdvtools spatial conversion against each SpatialData dataset "
            "in a directory and write a Markdown plus JSON report."
        )
    )
    parser.add_argument(
        "datasets_root",
        type=Path,
        help="Directory containing one SpatialData dataset per child entry.",
    )
    parser.add_argument(
        "report_dir",
        type=Path,
        help="Directory where the report, logs, and optional failure artifacts will be written.",
    )
    parser.add_argument(
        "--python-executable",
        type=Path,
        default=Path(sys.executable),
        help="Python executable used to invoke the conversion module.",
    )
    parser.add_argument(
        "--dataset",
        action="append",
        default=[],
        help="Specific dataset directory name to include. Repeat to include multiple.",
    )
    parser.add_argument(
        "--point-transform",
        choices=POINT_TRANSFORM_CHOICES,
        default="auto",
        help="Choose how table coordinates are transformed into image coordinates.",
    )
    parser.add_argument(
        "--output-geojson",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write transformed GeoJSON region files into the generated project.",
    )
    parser.add_argument(
        "--preserve-existing",
        action="store_true",
        help="Preserve existing project data in the output folder instead of recreating it.",
    )
    parser.add_argument(
        "--no-link",
        action="store_true",
        help="Copy the SpatialData inputs into outputs instead of linking them.",
    )
    parser.add_argument(
        "--keep-failures",
        action="store_true",
        help="Keep the per-dataset input/output temp directories for failed runs.",
    )
    parser.add_argument(
        "--density",
        action="store_true",
        help="Include density fields for gene expression in the default spatial view.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show detailed per-dataset conversion output, transform decisions, and merged summaries.",
    )
    return parser.parse_args()


def _discover_datasets(
    datasets_root: Path,
    selected_names: set[str],
    report_dir: Path | None = None,
) -> list[Path]:
    datasets = []
    resolved_report_dir = report_dir.resolve() if report_dir is not None else None
    for child in sorted(datasets_root.iterdir()):
        if child.name.startswith("."):
            continue
        if selected_names and child.name not in selected_names:
            continue
        if resolved_report_dir is not None and child.resolve() == resolved_report_dir:
            continue
        if child.is_dir():
            datasets.append(child)
    return datasets


def _ensure_dirs(report_dir: Path) -> tuple[Path, Path, Path]:
    logs_dir = report_dir / "logs"
    failures_dir = report_dir / "failures"
    cache_dir = report_dir / ".runtime-cache"
    for path in [report_dir, logs_dir, failures_dir, cache_dir]:
        path.mkdir(parents=True, exist_ok=True)
    return logs_dir, failures_dir, cache_dir


def _command_for_dataset(
    python_executable: Path,
    input_dir: Path,
    output_dir: Path,
    args: argparse.Namespace,
) -> list[str]:
    command = [
        str(python_executable),
        "-m",
        "mdvtools.spatial.conversion",
        str(input_dir),
        str(output_dir),
        "--point-transform",
        args.point_transform,
    ]
    if not args.no_link:
        command.append("--link")
    if args.output_geojson:
        command.append("--output_geojson")
    else:
        command.append("--no-output_geojson")
    if args.preserve_existing:
        command.append("--preserve-existing")
    if args.density:
        command.append("--density")
    if args.verbose:
        command.append("--verbose")
    return command


def _collect_warning_lines(stdout_text: str, stderr_text: str) -> list[str]:
    warnings: list[str] = []
    for line in [*stdout_text.splitlines(), *stderr_text.splitlines()]:
        if "warning" in line.lower():
            warnings.append(line.strip())
    return warnings


def _extract_traceback_excerpt(text: str) -> str | None:
    lines = text.splitlines()
    traceback_start = next(
        (index for index, line in enumerate(lines) if line.startswith("Traceback")),
        None,
    )
    if traceback_start is None:
        return None
    return "\n".join(lines[traceback_start:])


def _extract_error_details(
    returncode: int, stdout_text: str, stderr_text: str
) -> tuple[str | None, str | None, str | None, str | None]:
    if returncode == 0:
        return None, None, None, None

    traceback_excerpt = _extract_traceback_excerpt(stderr_text)
    if traceback_excerpt is None:
        traceback_excerpt = _extract_traceback_excerpt(stdout_text)
    if traceback_excerpt is not None:
        traceback_lines = [line.strip() for line in traceback_excerpt.splitlines() if line.strip()]
        last_traceback_line = traceback_lines[-1]
        error_type, separator, error_message = last_traceback_line.partition(":")
        if separator:
            signature = f"{error_type.strip()}: {error_message.strip()}"
            return (
                error_type.strip(),
                error_message.strip(),
                signature,
                traceback_excerpt,
            )
        return last_traceback_line, last_traceback_line, last_traceback_line, traceback_excerpt

    preferred_text = stderr_text.strip() or stdout_text.strip()
    non_empty_lines = [line.strip() for line in preferred_text.splitlines() if line.strip()]
    if not non_empty_lines:
        fallback = f"Process exited with code {returncode}"
        return None, fallback, fallback, None

    last_line = non_empty_lines[-1]
    error_type = None
    error_message = last_line
    signature = last_line
    if ":" in last_line:
        prefix, _, suffix = last_line.partition(":")
        if prefix.endswith("Error") or prefix.endswith("Exception"):
            error_type = prefix.strip()
            error_message = suffix.strip()
            signature = f"{error_type}: {error_message}"
    return error_type, error_message, signature, None


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _truncate_text(text: str, max_chars: int = 6000) -> str:
    if len(text) <= max_chars:
        return text
    return text[:max_chars].rstrip() + "\n...<traceback truncated>..."


def _run_dataset(
    dataset_path: Path,
    args: argparse.Namespace,
    logs_dir: Path,
    failures_dir: Path,
    cache_dir: Path,
) -> ConversionResult:
    dataset_name = dataset_path.name
    stdout_path = logs_dir / f"{dataset_name}.stdout.log"
    stderr_path = logs_dir / f"{dataset_name}.stderr.log"
    temp_root_path = Path(tempfile.mkdtemp(prefix=f"mdv-sdata-{dataset_name}-"))
    kept_failure_root = False
    try:
        input_dir = temp_root_path / "input"
        output_dir = temp_root_path / "output"
        input_dir.mkdir()
        output_dir.mkdir()
        os.symlink(dataset_path.resolve(), input_dir / dataset_name, target_is_directory=True)

        env = os.environ.copy()
        env["MPLCONFIGDIR"] = str(cache_dir / "matplotlib")
        env["XDG_CACHE_HOME"] = str(cache_dir / "xdg-cache")
        Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
        Path(env["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

        command = _command_for_dataset(
            args.python_executable, input_dir, output_dir, args
        )

        start = datetime.now(tz=UTC)
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )
        end = datetime.now(tz=UTC)

        stdout_text = completed.stdout
        stderr_text = completed.stderr
        _write_text(stdout_path, stdout_text)
        _write_text(stderr_path, stderr_text)

        if completed.returncode != 0 and args.keep_failures:
            kept_root = failures_dir / dataset_name
            if kept_root.exists():
                shutil.rmtree(kept_root)
            shutil.move(str(temp_root_path), kept_root)
            kept_failure_root = True
    finally:
        if temp_root_path.exists() and not kept_failure_root:
            shutil.rmtree(temp_root_path)

    warning_lines = _collect_warning_lines(stdout_text, stderr_text)
    error_type, error_message, error_signature, traceback_excerpt = _extract_error_details(
        completed.returncode,
        stdout_text,
        stderr_text,
    )

    return ConversionResult(
        dataset_name=dataset_name,
        dataset_path=str(dataset_path),
        status="success" if completed.returncode == 0 else "failed",
        returncode=completed.returncode,
        duration_seconds=round((end - start).total_seconds(), 2),
        command=command,
        stdout_log=str(stdout_path),
        stderr_log=str(stderr_path),
        warning_count=len(warning_lines),
        warning_samples=warning_lines[:5],
        error_type=error_type,
        error_message=error_message,
        error_signature=error_signature,
        traceback_excerpt=traceback_excerpt,
    )


def _format_markdown_report(
    datasets_root: Path,
    report_dir: Path,
    args: argparse.Namespace,
    results: list[ConversionResult],
) -> str:
    generated_at = datetime.now(tz=UTC).isoformat()
    success_count = sum(result.status == "success" for result in results)
    failure_count = len(results) - success_count
    failure_signatures = [
        result.error_signature
        for result in results
        if result.status == "failed" and result.error_signature is not None
    ]
    grouped_failures = Counter(failure_signatures)

    lines = [
        "# SpatialData conversion report",
        "",
        f"- Generated at: `{generated_at}`",
        f"- Datasets root: `{datasets_root}`",
        f"- Report directory: `{report_dir}`",
        f"- Python executable: `{args.python_executable}`",
        f"- Point transform: `{args.point_transform}`",
        f"- Link mode: `{not args.no_link}`",
        f"- Output geojson: `{args.output_geojson}`",
        f"- Datasets attempted: `{len(results)}`",
        f"- Successes: `{success_count}`",
        f"- Failures: `{failure_count}`",
        "",
        "## Summary",
        "",
        "| Dataset | Status | Seconds | Warnings | Error |",
        "| --- | --- | ---: | ---: | --- |",
    ]

    for result in results:
        error_summary = result.error_signature or ""
        lines.append(
            f"| `{result.dataset_name}` | `{result.status}` | `{result.duration_seconds}` | "
            f"`{result.warning_count}` | {error_summary} |"
        )

    if grouped_failures:
        lines.extend(
            [
                "",
                "## Failure groups",
                "",
            ]
        )
        for signature, count in grouped_failures.most_common():
            lines.append(f"- `{count}x` {signature}")

    lines.extend(["", "## Details", ""])
    for result in results:
        lines.extend(
            [
                f"### {result.dataset_name}",
                "",
                f"- Status: `{result.status}`",
                f"- Dataset path: `{result.dataset_path}`",
                f"- Duration seconds: `{result.duration_seconds}`",
                f"- Return code: `{result.returncode}`",
                f"- Command: `{shlex_join(result.command)}`",
                f"- Stdout log: `{result.stdout_log}`",
                f"- Stderr log: `{result.stderr_log}`",
            ]
        )
        if result.warning_samples:
            lines.append("- Warning samples:")
            for warning in result.warning_samples:
                lines.append(f"  - `{warning}`")
        if result.error_signature is not None:
            lines.append(f"- Error: `{result.error_signature}`")
        if result.traceback_excerpt is not None:
            lines.extend(
                [
                    "",
                    "```text",
                    _truncate_text(result.traceback_excerpt),
                    "```",
                ]
            )
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def shlex_join(parts: list[str]) -> str:
    return " ".join(_quote_shell_part(part) for part in parts)


def _quote_shell_part(part: str) -> str:
    if part == "":
        return "''"
    safe_characters = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-._/:=")
    if all(character in safe_characters for character in part):
        return part
    return "'" + part.replace("'", "'\"'\"'") + "'"


def main() -> int:
    args = _parse_args()
    datasets_root = args.datasets_root.expanduser().resolve()
    report_dir = args.report_dir.expanduser().resolve()
    selected_names = set(args.dataset)

    if not datasets_root.exists():
        raise FileNotFoundError(f"Datasets root does not exist: {datasets_root}")
    if not datasets_root.is_dir():
        raise NotADirectoryError(f"Datasets root is not a directory: {datasets_root}")

    logs_dir, failures_dir, cache_dir = _ensure_dirs(report_dir)
    datasets = _discover_datasets(datasets_root, selected_names, report_dir)
    if not datasets:
        raise ValueError(f"No datasets found under {datasets_root}")

    results: list[ConversionResult] = []
    for dataset in datasets:
        print(f"[run] {dataset.name}", flush=True)
        try:
            results.append(_run_dataset(dataset, args, logs_dir, failures_dir, cache_dir))
        except Exception as error:
            stdout_path = logs_dir / f"{dataset.name}.stdout.log"
            stderr_path = logs_dir / f"{dataset.name}.stderr.log"
            _write_text(stdout_path, "")
            _write_text(stderr_path, traceback.format_exc())
            results.append(
                ConversionResult(
                    dataset_name=dataset.name,
                    dataset_path=str(dataset),
                    status="failed",
                    returncode=-1,
                    duration_seconds=0.0,
                    command=[],
                    stdout_log=str(stdout_path),
                    stderr_log=str(stderr_path),
                    warning_count=0,
                    warning_samples=[],
                    error_type=type(error).__name__,
                    error_message=str(error),
                    error_signature=f"{type(error).__name__}: {error}",
                    traceback_excerpt=traceback.format_exc(),
                )
            )

    report_json_path = report_dir / "report.json"
    report_md_path = report_dir / "report.md"
    _write_text(
        report_json_path,
        json.dumps([asdict(result) for result in results], indent=2),
    )
    _write_text(
        report_md_path,
        _format_markdown_report(datasets_root, report_dir, args, results),
    )

    success_count = sum(result.status == "success" for result in results)
    failure_count = len(results) - success_count
    print(f"Report written to {report_md_path}")
    print(f"JSON written to {report_json_path}")
    print(f"Successes: {success_count}, failures: {failure_count}")
    return 0 if failure_count == 0 else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        traceback.print_exc()
        raise
