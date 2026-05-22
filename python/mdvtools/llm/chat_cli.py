import argparse
import csv
import json
import os
import re
import sys
import time
import traceback
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional
from urllib.parse import urlparse

from mdvtools.llm.chat_protocol import AskQuestionResult, ProjectChat
from mdvtools.mdvproject import MDVProject

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a one-shot Chat MDV prompt against an MDV project."
    )
    parser.add_argument("--project", required=True, help="Path to an MDV project folder.")
    parser.add_argument("--prompt", default=None, help="Natural-language prompt to run once.")
    parser.add_argument("--prompt-file", default=None, help="Path to newline-delimited prompts for batch mode.")
    parser.add_argument(
        "--view-name",
        default=None,
        help="Optional explicit name for the generated view.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Optional directory to store generated code and debug artifacts.",
    )
    parser.add_argument(
        "--csv-log",
        default=None,
        help="Path to CSV benchmark log output (required with --prompt-file).",
    )
    parser.add_argument(
        "--base-url",
        default=None,
        help="Base URL used to build full view URLs in batch CSV (required with --prompt-file).",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        dest="json_output",
        help="Emit machine-readable JSON output.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print extra diagnostic context in human-readable mode.",
    )
    return parser


class _CliSocketIOShim:
    def emit(self, _event: str, _payload: Any, namespace: Optional[str] = None, to: Optional[str] = None) -> None:
        del namespace, to
        return


def _extract_generated_code(code_blob: Optional[str]) -> Optional[str]:
    if not code_blob:
        return None
    match = re.search(r"```python\s*\n(.*?)```", code_blob, flags=re.DOTALL)
    if match:
        return match.group(1).strip()
    return code_blob.strip()


def _parse_block_timings(text: str) -> dict[str, float]:
    timings: dict[str, float] = {}
    pattern = re.compile(r"Block '([^']+)' took ([0-9]+(?:\.[0-9]+)?) seconds")
    for match in pattern.finditer(text):
        raw_name = match.group(1).strip().lower()
        safe_name = re.sub(r"[^a-z0-9]+", "_", raw_name).strip("_")
        if not safe_name:
            continue
        key = f"block_{safe_name}_s"
        timings[key] = float(match.group(2))
    return timings


def _count_charts_in_view(view: Any) -> int:
    """Count chart configs under MDV view[\"initialCharts\"]."""
    if not view or not isinstance(view, dict):
        return 0
    charts_by_ds = view.get("initialCharts")
    if not isinstance(charts_by_ds, dict):
        return 0
    n = 0
    for _ds, charts in charts_by_ds.items():
        if isinstance(charts, list):
            n += len(charts)
    return n


def _slugify(value: str) -> str:
    slug = re.sub(r"[^a-zA-Z0-9]+", "-", value).strip("-").lower()
    return slug[:40] or "prompt"


def _project_id_from_path(project_path: str) -> str:
    return Path(project_path).name


def _build_view_url(base_url: str, project_path: str, view_name: Optional[str]) -> str:
    if not view_name:
        return ""
    parsed = urlparse(base_url)
    if not parsed.scheme or not parsed.netloc:
        raise ValueError("--base-url must include scheme and host, e.g. http://localhost:5050")
    base = base_url.rstrip("/")
    return f"{base}/project/{_project_id_from_path(project_path)}/view/{view_name}"


def _normalize_result(
    *,
    project_path: str,
    views_file: str,
    ask_result: AskQuestionResult,
    debug_output_dir: Optional[str],
) -> dict[str, Any]:
    success = not ask_result.get("error", True)
    return {
        "success": success,
        "project_path": project_path,
        "views_file": views_file,
        "view_name": ask_result.get("view_name"),
        "message": ask_result.get("message", ""),
        "debug_output_dir": debug_output_dir,
    }


def _write_debug_artifacts(
    *,
    output_dir: str,
    ask_result: AskQuestionResult,
    normalized_result: dict[str, Any],
    view_snapshot: Any,
) -> str:
    out = Path(output_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    generated_code = _extract_generated_code(ask_result.get("code"))
    code_path = out / "generated_code.py"
    if generated_code:
        code_path.write_text(generated_code + "\n", encoding="utf-8")
    else:
        code_path.write_text("# No generated code available.\n", encoding="utf-8")

    result_payload = {
        "result": normalized_result,
        "raw_ask_result": ask_result,
    }
    (out / "result.json").write_text(
        json.dumps(result_payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    view_payload = {
        "view_name": normalized_result.get("view_name"),
        "view": view_snapshot,
    }
    (out / "view_snapshot.json").write_text(
        json.dumps(view_payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return str(out)


def run_chat_once(
    *,
    project_path: str,
    prompt: str,
    output_dir: Optional[str] = None,
    view_name: Optional[str] = None,
) -> tuple[dict[str, Any], int]:
    abs_project = str(Path(project_path).expanduser().resolve())
    views_file = str(Path(abs_project) / "views.json")
    debug_output_dir: Optional[str] = None
    timing_capture_start = time.perf_counter()

    try:
        captured_output = ""
        if not os.path.isdir(abs_project):
            raise FileNotFoundError(f"Project folder not found: {abs_project}")

        project = MDVProject(abs_project)
        if not project.datasources:
            raise ValueError("Project has no datasources; cannot run Chat MDV.")

        chat = ProjectChat(project)
        if getattr(chat, "init_error", False):
            raise RuntimeError(getattr(chat, "error_message", "Chat initialization failed."))

        from mdvtools import websocket as mdv_websocket
        if mdv_websocket.socketio is None:
            mdv_websocket.socketio = _CliSocketIOShim()  # type: ignore[assignment]

        def _handle_error(_error: Any, *, extra_metadata: Optional[dict] = None) -> None:
            del extra_metadata

        chat_request = {
            "message": prompt,
            "id": f"chat-cli-{uuid.uuid4().hex[:8]}",
            "conversation_id": f"chat-cli-{uuid.uuid4().hex[:8]}",
            "room": "chat-cli",
            "handle_error": _handle_error,
        }
        from io import StringIO
        from contextlib import redirect_stdout, redirect_stderr
        out_stream = StringIO()
        err_stream = StringIO()
        with redirect_stdout(out_stream), redirect_stderr(err_stream):
            ask_result = chat.ask_question(chat_request)
        captured_output = out_stream.getvalue() + "\n" + err_stream.getvalue()

        final_view_name = ask_result.get("view_name")
        if view_name:
            requested_view_name = view_name.strip()
            if not requested_view_name:
                raise ValueError("--view-name cannot be empty.")
            if final_view_name and requested_view_name != final_view_name:
                generated_view = project.get_view(final_view_name)
                if generated_view is not None:
                    project.set_view(requested_view_name, generated_view)
                    project.set_view(final_view_name, None)
                    final_view_name = requested_view_name
        ask_result["view_name"] = final_view_name
        view_snapshot = project.get_view(final_view_name) if final_view_name else None

        normalized = _normalize_result(
            project_path=abs_project,
            views_file=views_file,
            ask_result=ask_result,
            debug_output_dir=None,
        )
        normalized["block_timings"] = _parse_block_timings(captured_output)
        normalized["captured_output"] = captured_output
        normalized["duration_seconds"] = round(time.perf_counter() - timing_capture_start, 6)
        normalized["view_snapshot_present"] = view_snapshot is not None
        normalized["chart_count"] = _count_charts_in_view(view_snapshot)
        _ver = ask_result.get("verification")
        normalized["verification"] = _ver if isinstance(_ver, str) else ""
        normalized["needs_refresh"] = bool(ask_result.get("needs_refresh", False))

        if output_dir:
            debug_output_dir = _write_debug_artifacts(
                output_dir=output_dir,
                ask_result=ask_result,
                normalized_result=normalized,
                view_snapshot=view_snapshot,
            )
            normalized["debug_output_dir"] = debug_output_dir

        exit_code = 0 if normalized["success"] else 1
        return normalized, exit_code
    except Exception as exc:
        failure = {
            "success": False,
            "project_path": abs_project,
            "views_file": views_file,
            "view_name": None,
            "message": str(exc),
            "debug_output_dir": None,
            "block_timings": {},
            "captured_output": "",
            "duration_seconds": round(time.perf_counter() - timing_capture_start, 6),
            "view_snapshot_present": False,
            "chart_count": 0,
            "verification": "",
            "needs_refresh": False,
        }
        if output_dir:
            try:
                out = Path(output_dir).expanduser().resolve()
                out.mkdir(parents=True, exist_ok=True)
                failure["debug_output_dir"] = str(out)
                (out / "result.json").write_text(
                    json.dumps(
                        {
                            "result": failure,
                            "traceback": traceback.format_exc(),
                        },
                        indent=2,
                        sort_keys=True,
                    ),
                    encoding="utf-8",
                )
            except Exception:
                pass
        return failure, 2


def _validate_mode_args(args: argparse.Namespace) -> None:
    if bool(args.prompt) == bool(args.prompt_file):
        raise ValueError("Exactly one of --prompt or --prompt-file is required.")
    if args.prompt_file and not args.csv_log:
        raise ValueError("--csv-log is required when using --prompt-file.")
    if args.prompt_file and not args.base_url:
        raise ValueError("--base-url is required when using --prompt-file.")


def _load_prompts(prompt_file: str) -> list[str]:
    lines = Path(prompt_file).read_text(encoding="utf-8").splitlines()
    prompts = [line.strip() for line in lines if line.strip()]
    if not prompts:
        raise ValueError("Prompt file is empty.")
    return prompts


def _write_csv_rows(csv_path: str, rows: list[dict[str, Any]]) -> None:
    fixed_fields = [
        "row_index",
        "prompt",
        "status",
        "exit_code",
        "error_message",
        "started_at",
        "finished_at",
        "duration_seconds",
        "project_path",
        "view_name",
        "view_url",
        "code_path",
        "result_json_path",
        "view_snapshot_path",
    ]
    timing_keys: set[str] = set()
    for row in rows:
        for key in row:
            if key.startswith("block_") and key.endswith("_s"):
                timing_keys.add(key)
    fieldnames = fixed_fields + sorted(timing_keys)
    out = Path(csv_path).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def run_batch_prompts(
    *,
    project_path: str,
    prompt_file: str,
    csv_log: str,
    base_url: str,
    output_dir: Optional[str] = None,
) -> tuple[list[dict[str, Any]], int]:
    prompts = _load_prompts(prompt_file)
    rows: list[dict[str, Any]] = []
    failures = 0
    base_out = Path(output_dir).expanduser().resolve() if output_dir else None
    if base_out:
        base_out.mkdir(parents=True, exist_ok=True)

    for idx, prompt in enumerate(prompts, start=1):
        started = datetime.now(timezone.utc)
        run_out = None
        if base_out:
            run_out = base_out / f"run_{idx:04d}_{_slugify(prompt)}"
        result, exit_code = run_chat_once(
            project_path=project_path,
            prompt=prompt,
            output_dir=str(run_out) if run_out else None,
            view_name=None,
        )
        finished = datetime.now(timezone.utc)
        if exit_code != 0:
            failures += 1
        debug_dir = result.get("debug_output_dir")
        code_path = str(Path(debug_dir) / "generated_code.py") if debug_dir else ""
        result_json_path = str(Path(debug_dir) / "result.json") if debug_dir else ""
        view_snapshot_path = str(Path(debug_dir) / "view_snapshot.json") if debug_dir else ""

        row: dict[str, Any] = {
            "row_index": idx,
            "prompt": prompt,
            "status": "success" if result.get("success") else "failure",
            "exit_code": exit_code,
            "error_message": "" if result.get("success") else result.get("message", ""),
            "started_at": started.isoformat(),
            "finished_at": finished.isoformat(),
            "duration_seconds": result.get("duration_seconds", 0.0),
            "project_path": result.get("project_path", project_path),
            "view_name": result.get("view_name", ""),
            "view_url": _build_view_url(base_url, project_path, result.get("view_name")),
            "code_path": code_path,
            "result_json_path": result_json_path,
            "view_snapshot_path": view_snapshot_path,
        }
        for k, v in (result.get("block_timings") or {}).items():
            row[k] = v
        rows.append(row)

    _write_csv_rows(csv_log, rows)
    return rows, (0 if failures == 0 else 1)


def _print_human_result(result: dict[str, Any], *, verbose: bool = False) -> None:
    status = "SUCCESS" if result["success"] else "FAILED"
    print(status)
    print(f"Message: {result['message']}")
    print(f"Project: {result['project_path']}")
    print(f"Views file: {result['views_file']}")
    print(f"View: {result['view_name']}")
    if result.get("debug_output_dir"):
        print(f"Debug output dir: {result['debug_output_dir']}")
    elif verbose:
        print("Debug output dir: None")


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    try:
        args = parser.parse_args(argv)
        _validate_mode_args(args)
        if args.prompt_file:
            rows, exit_code = run_batch_prompts(
                project_path=args.project,
                prompt_file=args.prompt_file,
                csv_log=args.csv_log,
                base_url=args.base_url,
                output_dir=args.output_dir,
            )
            summary = {
                "batch_count": len(rows),
                "success_count": sum(1 for r in rows if r.get("status") == "success"),
                "failure_count": sum(1 for r in rows if r.get("status") == "failure"),
                "csv_log": str(Path(args.csv_log).expanduser().resolve()),
            }
            if args.json_output:
                print(json.dumps(summary, indent=2, sort_keys=True))
            else:
                print(f"Batch completed: {summary['success_count']} success, {summary['failure_count']} failure")
                print(f"CSV log: {summary['csv_log']}")
            return exit_code

        result, exit_code = run_chat_once(
            project_path=args.project,
            prompt=args.prompt,
            output_dir=args.output_dir,
            view_name=args.view_name,
        )

        if args.json_output:
            print(json.dumps(result, indent=2, sort_keys=True))
        else:
            _print_human_result(result, verbose=args.verbose)

        return exit_code
    except ValueError as exc:
        print(f"Argument error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
