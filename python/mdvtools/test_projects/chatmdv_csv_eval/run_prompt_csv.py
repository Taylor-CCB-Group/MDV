#!/usr/bin/env python3
"""
Run Chat MDV prompts from a CSV (one row per benchmark), using each row's project path.

Equivalent pipeline to ``python -m mdvtools.llm.chat_cli`` with ``--prompt`` /
``--project`` via ``chat_cli.run_chat_once``.

Setup
-----

From the ``python`` directory where package ``mdvtools`` resolves (repository
``python/``):

  python mdvtools/test_projects/chatmdv_csv_eval/run_prompt_csv.py \
    --csv mdvtools/test_projects/chatmdv_csv_eval/prompts_pbmc3k.csv

Each row ``Path`` must be an existing MDV project folder (with ``datasources``, etc.).
Example in this workspace: ``/app/mdv/pbmc3k_chat`` (used in ``prompts_pbmc3k.csv``).

Environment: same as Chat MDV (e.g. ``OPENAI_API_KEY`` in ``.env`` if used).

Default report directory for outputs: ``$TMPDIR/chatmdv_csv_eval`` (on Linux usually
under ``/tmp``). Override with ``--output-dir /path``, env ``CHATMDV_CSV_EVAL_OUTPUT_DIR``,
or an explicit ``--results`` / ``--failure-log`` / ``--details-log`` path.

Input CSV columns (flexible names)
----------------------------------

- **Path** (or ``path``, ``project_path``): MDV project directory.
- **Question** (or ``Question for evaluation``, ``prompt``): natural-language prompt.
- Optional: **Dataset**, **Complexity score** (or ``complexity``)—copied through and
  included in detail/failure records.

Outputs
-------

- **Results CSV** (default ``results_<stem>_<utc_timestamp>.csv``): original columns plus
  compact per-row metrics only:

  - ``exec_success``, ``exit_code``, ``duration_seconds``
  - ``view_name``, ``view_snapshot_present``, ``chart_count``, ``view_has_charts``,
    ``needs_refresh``
  - ``details_jsonl_path`` — path to the sidecar details file for this run

  Verbose fields are **not** in the results CSV (see details JSONL below).

- **Details JSONL** (``details_<stem>_<timestamp>.jsonl``): one JSON object per processed
  row with ``failure_reason``, ``captured_output_excerpt``, ``captured_output``,
  ``verification``, ``debug_output_dir``, and row identity fields.

- **Benchmark summary JSON** (``benchmark_summary_<stem>_<timestamp>.json``): run metadata
  (``run_id``, ``run_completed_at_utc``, ``run_wall_seconds_total``, ``run_cli_limit``,
  ``run_artifacts_base_dir``), artifact paths, ``totals`` (including ``worker_sigkill_rows``
  and ``row_isolation_subprocess``), and ``by_complexity``.

- **Failures JSONL** (``failures_<stem>_<timestamp>.jsonl``): written only if at least one
  row failed; same schema as a failed row in the details JSONL (quick filter).

Row isolation (OOM recovery)
------------------------------

By default each benchmark row runs in a **child Python process**. If the OS OOM killer
stops that child with **SIGKILL**, the parent records ``exit_code`` **3** with a
``failure_reason`` explaining likely OOM and **continues** with the next row. Use
``--in-process-rows`` to run everything in one process (faster imports, but a single OOM
kill ends the whole CSV run). Worker subprocess failures before returning JSON use
``exit_code`` **4**.

Use ``--artifacts-dir`` to persist per-row ``generated_code.py`` / ``result.json``
(see ``mdvtools.llm.chat_cli``).

Per-row ``exit_code`` is ``chat_cli.run_chat_once`` semantics (**0** success, **1**
logical failure, **2** exception) when the worker returns normally; **3** worker killed
by signal (SIGKILL is the usual Linux OOM case), **4** worker handshake/bootstrap failure,
**5** worker subprocess timed out (see ``--row-timeout`` / ``ROW_TIMEOUT_SECONDS``).
Process exit: **0** if all rows succeeded, **1** if any row failed, **2** only for bad
paths or empty CSV header.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

_PY_ROOT = Path(__file__).resolve().parents[3]
_SCRIPT_PATH = Path(__file__).resolve()
if str(_PY_ROOT) not in sys.path:
    sys.path.insert(0, str(_PY_ROOT))

# SIGKILL is not always exposed on all platforms (e.g. some Windows builds).
_SIGKILL = getattr(signal, "SIGKILL", 9)


EXCERPT_MAX_LEN = 4000

_DEFAULT_ROW_TIMEOUT_RAW = os.environ.get("CHATMDV_CSV_EVAL_ROW_TIMEOUT_SECONDS", "3600").strip()
try:
    ROW_TIMEOUT_SECONDS: float = float(_DEFAULT_ROW_TIMEOUT_RAW)
except ValueError:
    ROW_TIMEOUT_SECONDS = 3600.0

SLIM_RESULT_FIELDS = [
    "exec_success",
    "exit_code",
    "duration_seconds",
    "view_name",
    "view_snapshot_present",
    "chart_count",
    "view_has_charts",
    "needs_refresh",
    "details_jsonl_path",
]


def _default_report_root() -> Path:
    raw = os.environ.get("CHATMDV_CSV_EVAL_OUTPUT_DIR", "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    return (Path(tempfile.gettempdir()) / "chatmdv_csv_eval").resolve()


def _norm_key(h: str) -> str:
    return h.strip().lower().replace(" ", "_")


def _column_aliases() -> dict[str, list[str]]:
    return {
        "path": ["path", "project_path"],
        "question": ["question_for_evaluation", "question", "prompt"],
        "dataset": ["dataset"],
        "complexity": ["complexity_score", "complexity"],
    }


def _resolve_row_field(row: dict[str, str], logical: str) -> str | None:
    aliases = set(_column_aliases().get(logical, []))
    for hk, val in row.items():
        nk = _norm_key(hk)
        if nk in aliases:
            v = (val or "").strip()
            if v:
                return v
    return None


def _excerpt_tail(text: str, max_len: int = EXCERPT_MAX_LEN) -> str:
    if not text:
        return ""
    if len(text) <= max_len:
        return text
    return f"<{len(text) - max_len} chars omitted>\n{text[-max_len:]}"


def _row_detail_payload(
    *,
    row_index: int,
    path: str,
    dataset: str,
    question: str,
    complexity_score: str,
    exec_success: bool,
    exit_code: int,
    duration_seconds: float,
    failure_reason: str,
    captured_output: str,
    captured_output_excerpt: str,
    view_name: str,
    view_snapshot_present: bool,
    chart_count: int,
    view_has_charts: bool,
    verification: str,
    needs_refresh: bool,
    debug_output_dir: str,
) -> dict[str, object]:
    return {
        "row_index": row_index,
        "path": path,
        "dataset": dataset,
        "question": question,
        "complexity_score": complexity_score,
        "exec_success": exec_success,
        "exit_code": exit_code,
        "duration_seconds": duration_seconds,
        "failure_reason": failure_reason,
        "captured_output_excerpt": captured_output_excerpt,
        "captured_output": captured_output,
        "view_name": view_name,
        "view_snapshot_present": view_snapshot_present,
        "chart_count": chart_count,
        "view_has_charts": view_has_charts,
        "verification": verification,
        "needs_refresh": needs_refresh,
        "debug_output_dir": debug_output_dir,
    }


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--csv", "-c", required=True, type=Path, help="Benchmark CSV path.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for results, benchmark summary JSON, details/failure JSONL (default: "
        "$TMPDIR/chatmdv_csv_eval, or CHATMDV_CSV_EVAL_OUTPUT_DIR). Ignored if --results is set.",
    )
    parser.add_argument(
        "--results",
        "-o",
        type=Path,
        default=None,
        help="Results CSV output path (default: under output-dir/report root with timestamp).",
    )
    parser.add_argument(
        "--details-log",
        type=Path,
        default=None,
        help="JSONL path for per-row verbose details (default: details_* next to results CSV).",
    )
    parser.add_argument(
        "--failure-log",
        type=Path,
        default=None,
        help="JSONL path for failures (default: same directory as results, failures_* prefix).",
    )
    parser.add_argument(
        "--artifacts-dir",
        type=Path,
        default=None,
        help="Base directory for per-row debug dirs (passed to chat_cli.run_chat_once).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process at most N data rows after the header (for smoke tests).",
    )
    parser.add_argument(
        "--in-process-rows",
        action="store_true",
        help=(
            "Run each CSV row in the main Python process (lower overhead). "
            "If the OS OOM-kills the process, the entire CSV run stops; default is one "
            "worker subprocess per row so SIGKILL is recorded and the run continues."
        ),
    )
    parser.add_argument(
        "--row-timeout",
        type=float,
        default=None,
        metavar="SECONDS",
        help=(
            "Per-row worker subprocess wall-clock timeout (default: ROW_TIMEOUT_SECONDS "
            f"env or {ROW_TIMEOUT_SECONDS:g}s). Ignored with --in-process-rows."
        ),
    )
    return parser.parse_args()


def _find_header(fieldnames: list[str], logical: str) -> str | None:
    aliases = {a.replace(" ", "_") for a in _column_aliases().get(logical, [])}
    for fn in fieldnames:
        if _norm_key(fn) in aliases:
            return fn
    return None


def _slug(text: str, max_len: int = 40) -> str:
    s = re.sub(r"[^a-zA-Z0-9]+", "-", text).strip("-").lower()
    return (s[:max_len] or "prompt").rstrip("-")


def _complexity_sort_key(k: str) -> tuple[int, str]:
    try:
        return (0, f"{int(k):010d}")
    except ValueError:
        return (1, k)


def _synthetic_worker_failure(
    *,
    project_path: str,
    message: str,
    stderr_excerpt: str = "",
) -> dict[str, object]:
    cap = stderr_excerpt.strip()
    return {
        "success": False,
        "project_path": project_path,
        "views_file": str(Path(project_path).expanduser().resolve() / "views.json"),
        "view_name": None,
        "message": message,
        "debug_output_dir": None,
        "block_timings": {},
        "captured_output": cap,
        "duration_seconds": 0.0,
        "view_snapshot_present": False,
        "chart_count": 0,
        "verification": "",
        "needs_refresh": False,
    }


def _signal_label(sig_num: int) -> str:
    if hasattr(signal, "strsignal"):
        try:
            label = signal.strsignal(sig_num)
            if label:
                return label
        except Exception:
            pass
    return f"signal {sig_num}"


def _csv_eval_worker_main() -> int:
    """stdin: one JSON object with project_path, prompt, optional output_dir. stdout: one JSON line."""
    try:
        raw_in = sys.stdin.read()
        payload = json.loads(raw_in)
    except Exception as exc:
        sys.stdout.write(
            json.dumps({"worker_bootstrap_error": str(exc), "traceback": traceback.format_exc()})
            + "\n"
        )
        sys.stdout.flush()
        return 1

    project_path = str(payload.get("project_path") or "")
    prompt = str(payload.get("prompt") or "")
    output_dir = payload.get("output_dir")
    out_dir: str | None = str(output_dir) if output_dir else None

    from mdvtools.llm.chat_cli import run_chat_once

    try:
        result, exit_code = run_chat_once(
            project_path=project_path,
            prompt=prompt,
            output_dir=out_dir,
            view_name=None,
        )
    except Exception as exc:
        sys.stdout.write(
            json.dumps(
                {
                    "worker_error": str(exc),
                    "traceback": traceback.format_exc(),
                    "project_path": project_path,
                }
            )
            + "\n"
        )
        sys.stdout.flush()
        return 1

    sys.stdout.write(json.dumps({"result": result, "exit_code": exit_code}, ensure_ascii=False) + "\n")
    sys.stdout.flush()
    return 0


def _run_chat_row_subprocess(
    *,
    project_path: str,
    prompt: str,
    output_dir: str | None,
    timeout_seconds: float | None = None,
) -> tuple[dict[str, object], int, str]:
    """Run ``run_chat_once`` in a child process.

    Returns ``(result_dict, row_exit_code, termination_tag)`` where ``termination_tag`` is
    ``\"\"`` for a normal worker handshake, ``\"sigkill\"`` when the child was SIGKILL'd
    (common for Linux OOM), ``\"signal_other\"`` for another fatal signal, or
    ``\"timeout\"`` when the worker exceeded ``timeout_seconds``.
    """
    payload = {"project_path": project_path, "prompt": prompt, "output_dir": output_dir}
    cmd = [sys.executable, str(_SCRIPT_PATH), "--csv-eval-worker"]
    run_kwargs: dict[str, object] = {
        "input": json.dumps(payload, ensure_ascii=False).encode("utf-8"),
        "capture_output": True,
        "check": False,
        "cwd": str(_PY_ROOT),
    }
    if timeout_seconds is not None and timeout_seconds > 0:
        run_kwargs["timeout"] = timeout_seconds
    try:
        completed = subprocess.run(cmd, **run_kwargs)
    except subprocess.TimeoutExpired as exc:
        stderr_txt = (exc.stderr or b"").decode("utf-8", errors="replace")
        stdout_txt = (exc.stdout or b"").decode("utf-8", errors="replace")
        cap_parts = []
        if stdout_txt.strip():
            cap_parts.append("Worker stdout excerpt:\n" + _excerpt_tail(stdout_txt, EXCERPT_MAX_LEN))
        if stderr_txt.strip():
            cap_parts.append("Worker stderr excerpt:\n" + _excerpt_tail(stderr_txt, EXCERPT_MAX_LEN))
        cap = "\n".join(cap_parts)
        msg = (
            f"Worker subprocess timed out after {timeout_seconds:g}s; "
            "the CSV runner continues with the next row."
        )
        fail = _synthetic_worker_failure(
            project_path=project_path,
            message=msg,
            stderr_excerpt=cap or stderr_txt,
        )
        return fail, 5, "timeout"
    except OSError as exc:
        fail = _synthetic_worker_failure(
            project_path=project_path,
            message=f"Failed to spawn CSV eval worker subprocess: {exc}",
        )
        return fail, 4, ""

    rc = completed.returncode
    stderr_txt = (completed.stderr or b"").decode("utf-8", errors="replace")
    stdout_txt = (completed.stdout or b"").decode("utf-8", errors="replace").strip()

    if rc < 0:
        sig = -rc
        if sig == _SIGKILL:
            msg = (
                "Worker subprocess was killed with SIGKILL (signal 9). On Linux this "
                "commonly indicates the OOM killer stopped the child after memory pressure; "
                "the CSV runner continues with the next row."
            )
            tag = "sigkill"
        else:
            msg = f"Worker subprocess was killed ({_signal_label(sig)})."
            tag = "signal_other"
        if stderr_txt.strip():
            msg += " Worker stderr excerpt:\n" + _excerpt_tail(stderr_txt, EXCERPT_MAX_LEN)
        fail = _synthetic_worker_failure(project_path=project_path, message=msg, stderr_excerpt=stderr_txt)
        return fail, 3, tag

    def _parse_worker_stdout(text: str) -> dict[str, object] | None:
        if not text:
            return None
        for line in reversed(text.splitlines()):
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(obj, dict):
                return obj
        return None

    parsed = _parse_worker_stdout(stdout_txt)

    if rc != 0:
        hint = ""
        if isinstance(parsed, dict):
            if "worker_bootstrap_error" in parsed:
                hint = str(parsed.get("worker_bootstrap_error") or "")
            elif "worker_error" in parsed:
                hint = str(parsed.get("worker_error") or "")
        msg = f"Worker subprocess exited with code {rc} before returning a usable result."
        if hint:
            msg += f" Detail: {hint}"
        if stderr_txt.strip():
            msg += "\nWorker stderr excerpt:\n" + _excerpt_tail(stderr_txt, EXCERPT_MAX_LEN)
        fail = _synthetic_worker_failure(project_path=project_path, message=msg, stderr_excerpt=stderr_txt)
        return fail, 4, ""

    if not isinstance(parsed, dict) or "result" not in parsed or "exit_code" not in parsed:
        excerpt = _excerpt_tail(stdout_txt, 500)
        msg = (
            "Worker subprocess returned stdout that is not a valid result payload "
            f"(stdout length {len(stdout_txt)}; excerpt): {excerpt!r}"
        )
        fail = _synthetic_worker_failure(project_path=project_path, message=msg, stderr_excerpt=stderr_txt)
        return fail, 4, ""

    result_obj = parsed["result"]
    if not isinstance(result_obj, dict):
        fail = _synthetic_worker_failure(
            project_path=project_path,
            message="Worker JSON payload field 'result' is not an object.",
            stderr_excerpt=stderr_txt,
        )
        return fail, 4, ""

    try:
        row_exit = int(parsed["exit_code"])
    except (TypeError, ValueError):
        fail = _synthetic_worker_failure(
            project_path=project_path,
            message=f"Worker JSON payload has invalid exit_code: {parsed.get('exit_code')!r}",
            stderr_excerpt=stderr_txt,
        )
        return fail, 4, ""

    return result_obj, row_exit, ""


def main() -> int:
    args = _parse_args()
    isolate_rows = not bool(args.in_process_rows)
    row_timeout_seconds = (
        args.row_timeout if args.row_timeout is not None else ROW_TIMEOUT_SECONDS
    )

    csv_path = args.csv.expanduser().resolve()
    if not csv_path.is_file():
        print(f"Not found: {csv_path}", file=sys.stderr)
        return 2

    stem = csv_path.stem
    utc_now = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    if args.results is None:
        report_root = (
            args.output_dir.expanduser().resolve() if args.output_dir else _default_report_root()
        )
        report_root.mkdir(parents=True, exist_ok=True)
        results_path = report_root / f"results_{stem}_{utc_now}.csv"
    else:
        results_path = args.results.expanduser().resolve()
    results_path.parent.mkdir(parents=True, exist_ok=True)
    benchmark_summary_path = results_path.parent / f"benchmark_summary_{stem}_{utc_now}.json"

    details_log_path = args.details_log
    if details_log_path is None:
        details_log_path = results_path.parent / f"details_{stem}_{utc_now}.jsonl"
    else:
        details_log_path = details_log_path.expanduser().resolve()
    details_log_path.parent.mkdir(parents=True, exist_ok=True)
    details_log_path_str = str(details_log_path)

    failure_log_path = args.failure_log
    if failure_log_path is None:
        failure_log_path = results_path.parent / f"failures_{stem}_{utc_now}.jsonl"
    else:
        failure_log_path = failure_log_path.expanduser().resolve()
    failure_log_path.parent.mkdir(parents=True, exist_ok=True)

    base_artifacts = args.artifacts_dir
    if base_artifacts is not None:
        base_artifacts = base_artifacts.expanduser().resolve()
        base_artifacts.mkdir(parents=True, exist_ok=True)

    # utf-8-sig strips a leading UTF-8 BOM so the first column is "Path", not "\ufeffPath"
    with csv_path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            print("CSV has no header row.", file=sys.stderr)
            return 2
        raw_fieldnames = list(reader.fieldnames)
        rows_in = list(reader)

    if args.limit is not None:
        rows_in = rows_in[: args.limit]

    wall_start = time.perf_counter()

    result_rows: list[dict[str, object]] = []
    detail_records: list[dict[str, object]] = []
    complexity_pairs: dict[str, list[tuple[bool, bool]]] = defaultdict(list)
    failure_count = 0
    failure_payloads: list[dict[str, object]] = []
    worker_sigkill_rows = 0

    comp_header = _find_header(raw_fieldnames, "complexity")
    n_rows = len(rows_in)

    timeout_note = (
        f", row timeout={row_timeout_seconds:g}s"
        if isolate_rows and row_timeout_seconds > 0
        else ""
    )
    print(
        f"Processing {n_rows} row(s) from {csv_path.name} "
        f"(row isolation subprocess={isolate_rows}{timeout_note})",
        file=sys.stderr,
        flush=True,
    )

    for idx, row in enumerate(rows_in, start=1):
        path = _resolve_row_field(row, "path")
        question = _resolve_row_field(row, "question")
        dataset = _resolve_row_field(row, "dataset") or ""
        complexity_val = ((row.get(comp_header) or "").strip() if comp_header else "") or ""

        q_preview = (question or "")[:80] + ("..." if question and len(question) > 80 else "")
        print(
            f"[START] row={idx}/{n_rows} path={path or '(missing)'} q={q_preview!r}",
            file=sys.stderr,
            flush=True,
        )

        if not path or not question:
            msg = f"missing Path or Question (path={path!r}, question={question!r})"
            print(f"[FAIL] row={idx} {msg}", file=sys.stderr, flush=True)
            detail = _row_detail_payload(
                row_index=idx,
                path=path or "",
                dataset=dataset,
                question=question or "",
                complexity_score=complexity_val,
                exec_success=False,
                exit_code=2,
                duration_seconds=0.0,
                failure_reason=f"Row {idx}: {msg}",
                captured_output="",
                captured_output_excerpt="",
                view_name="",
                view_snapshot_present=False,
                chart_count=0,
                view_has_charts=False,
                verification="",
                needs_refresh=False,
                debug_output_dir="",
            )
            detail_records.append(detail)
            out_row = {
                **{k: row.get(k, "") for k in raw_fieldnames},
                "exec_success": False,
                "exit_code": 2,
                "duration_seconds": 0.0,
                "view_name": "",
                "view_snapshot_present": False,
                "chart_count": 0,
                "view_has_charts": False,
                "needs_refresh": False,
                "details_jsonl_path": details_log_path_str,
            }
            result_rows.append(out_row)
            failure_count += 1
            if complexity_val:
                complexity_pairs[complexity_val].append((False, False))
            failure_payloads.append(detail)
            continue

        run_out = None
        if base_artifacts is not None:
            run_out = base_artifacts / f"row_{idx:04d}_{_slug(question)}"

        out_dir_str = str(run_out) if run_out else None
        term_tag = ""
        if isolate_rows:
            row_timeout = row_timeout_seconds if row_timeout_seconds > 0 else None
            result, exit_code, term_tag = _run_chat_row_subprocess(
                project_path=path,
                prompt=question,
                output_dir=out_dir_str,
                timeout_seconds=row_timeout,
            )
            if term_tag == "timeout":
                print(
                    f"[TIMEOUT] row={idx}: worker subprocess exceeded {row_timeout_seconds:g}s; "
                    f"recorded as exit_code=5; continuing.",
                    file=sys.stderr,
                    flush=True,
                )
            elif term_tag == "sigkill":
                worker_sigkill_rows += 1
                print(
                    f"[OOM?] row={idx}: worker subprocess was SIGKILL'd (Linux OOM killer "
                    f"often uses SIGKILL). Recorded as exit_code=3; continuing.",
                    file=sys.stderr,
                    flush=True,
                )
            elif term_tag == "signal_other":
                print(
                    f"[WORKER] row={idx}: worker subprocess was killed by a signal other than "
                    f"SIGKILL; recorded as exit_code=3; continuing.",
                    file=sys.stderr,
                    flush=True,
                )
        else:
            from mdvtools.llm.chat_cli import run_chat_once as run_chat_once_fn

            result, exit_code = run_chat_once_fn(
                project_path=path,
                prompt=question,
                output_dir=out_dir_str,
                view_name=None,
            )
        success = bool(result.get("success"))
        message = str(result.get("message", ""))
        cap = str(result.get("captured_output", ""))
        duration = float(result.get("duration_seconds", 0.0))
        view_name = result.get("view_name") or ""
        debug_dir = result.get("debug_output_dir") or ""
        view_snapshot_present = bool(result.get("view_snapshot_present", False))
        chart_count = int(result.get("chart_count") or 0)
        _ver = result.get("verification")
        verification = _ver if isinstance(_ver, str) else ""
        needs_refresh = bool(result.get("needs_refresh", False))
        view_has_charts = success and chart_count > 0

        failure_reason = "" if success else message
        excerpt = "" if success else _excerpt_tail(cap)

        detail = _row_detail_payload(
            row_index=idx,
            path=path,
            dataset=dataset,
            question=question,
            complexity_score=complexity_val,
            exec_success=success,
            exit_code=exit_code,
            duration_seconds=duration,
            failure_reason=failure_reason,
            captured_output=cap,
            captured_output_excerpt=excerpt,
            view_name=view_name,
            view_snapshot_present=view_snapshot_present,
            chart_count=chart_count,
            view_has_charts=view_has_charts,
            verification=verification,
            needs_refresh=needs_refresh,
            debug_output_dir=str(debug_dir) if debug_dir else "",
        )
        detail_records.append(detail)

        out_row = {
            **{k: row.get(k, "") for k in raw_fieldnames},
            "exec_success": success,
            "exit_code": exit_code,
            "duration_seconds": duration,
            "view_name": view_name,
            "view_snapshot_present": view_snapshot_present,
            "chart_count": chart_count,
            "view_has_charts": view_has_charts,
            "needs_refresh": needs_refresh,
            "details_jsonl_path": details_log_path_str,
        }
        result_rows.append(out_row)

        if complexity_val:
            complexity_pairs[complexity_val].append((success, view_has_charts))

        if success:
            view_bit = f" view={view_name!r}" if view_name else ""
            charts_bit = f" charts={chart_count}" if view_has_charts else ""
            print(
                f"[OK] row={idx}/{n_rows} exit={exit_code} duration={duration:.1f}s"
                f"{view_bit}{charts_bit}",
                file=sys.stderr,
                flush=True,
            )
        else:
            failure_count += 1
            q_short = question[:120] + ("..." if len(question) > 120 else "")
            print(
                f"[FAIL] row={idx}/{n_rows} exit={exit_code} path={path} q={q_short!r}\n"
                f"       reason: {failure_reason}",
                file=sys.stderr,
                flush=True,
            )
            failure_payloads.append(detail)

    with details_log_path.open("w", encoding="utf-8") as dlog:
        for record in detail_records:
            dlog.write(json.dumps(record, ensure_ascii=False) + "\n")

    if failure_payloads:
        with failure_log_path.open("w", encoding="utf-8") as flog:
            for payload in failure_payloads:
                flog.write(json.dumps(payload, ensure_ascii=False) + "\n")

    n = len(result_rows)
    ok = n - failure_count
    rate = (ok / n * 100.0) if n else 0.0
    view_charts_ok = sum(1 for r in result_rows if r.get("view_has_charts"))
    view_rate = (view_charts_ok / n * 100.0) if n else 0.0

    by_complexity: dict[str, dict[str, int]] = {}
    for key in sorted(complexity_pairs.keys(), key=_complexity_sort_key):
        pairs = complexity_pairs[key]
        c_n = len(pairs)
        exec_ok = sum(1 for s, _ in pairs if s)
        vch_ok = sum(1 for _, v in pairs if v)
        by_complexity[str(key)] = {
            "rows": c_n,
            "exec_success_count": exec_ok,
            "view_has_charts_count": vch_ok,
        }

    completion_iso = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")
    wall_seconds = round(time.perf_counter() - wall_start, 6)

    summary_payload = {
        "run_id": utc_now,
        "run_completed_at_utc": completion_iso,
        "run_wall_seconds_total": wall_seconds,
        "run_cli_limit": args.limit,
        "run_artifacts_base_dir": str(base_artifacts) if base_artifacts else None,
        "input_csv": str(csv_path),
        "results_csv": str(results_path),
        "details_jsonl": str(details_log_path),
        "failure_log": str(failure_log_path) if failure_count > 0 else None,
        "benchmark_summary": str(benchmark_summary_path),
        "totals": {
            "rows": n,
            "exec_success_count": ok,
            "exec_success_rate": round(ok / n, 6) if n else 0.0,
            "view_has_charts_count": view_charts_ok,
            "view_has_charts_rate": round(view_charts_ok / n, 6) if n else 0.0,
            "worker_sigkill_rows": worker_sigkill_rows,
            "row_isolation_subprocess": isolate_rows,
        },
        "by_complexity": by_complexity,
    }

    out_fieldnames = list(raw_fieldnames)
    for f in SLIM_RESULT_FIELDS:
        if f not in out_fieldnames:
            out_fieldnames.append(f)

    with results_path.open("w", newline="", encoding="utf-8") as wf:
        w = csv.DictWriter(wf, fieldnames=out_fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in result_rows:
            w.writerow({k: r.get(k, "") for k in out_fieldnames})

    benchmark_summary_path.write_text(
        json.dumps(summary_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    report_dir = results_path.parent.resolve()
    print(f"Processed {n} row(s): {ok} success, {failure_count} failure ({rate:.1f}% success rate)")
    print(f"View with charts (structural): {view_charts_ok}/{n} ({view_rate:.1f}%)")
    if isolate_rows and worker_sigkill_rows:
        print(
            f"Worker SIGKILL rows (likely OOM; see per-row failure_reason / exit_code=3): "
            f"{worker_sigkill_rows}",
            file=sys.stderr,
        )

    print("Output directory (persisted artifacts):")
    print(f"  {report_dir}")
    print("Files written:")
    print(f"  results_csv               {results_path}")
    print(f"  details_jsonl             {details_log_path}")
    print(f"  benchmark_summary_json    {benchmark_summary_path}")
    if failure_count:
        print(f"  failures_jsonl            {failure_log_path}")
    else:
        print("  failures_jsonl            (not written — no failed rows)")
    if base_artifacts:
        print(f"  artifacts_dir (--artifacts-dir) {base_artifacts}")
    else:
        print("  artifacts_dir             (none — pass --artifacts-dir for per-row debug dirs)")
    print("")
    print(
        "During the run: [START]/[OK]/[FAIL] lines on stderr (flushed per row)."
        + (
            f' Verbose per-row fields: see details_jsonl ("{details_log_path.name}").'
            f' Failed-row payloads also in failures_jsonl ("{failure_log_path.name}").'
            if failure_count
            else f' Verbose per-row fields: see details_jsonl ("{details_log_path.name}").'
        )
    )

    if complexity_pairs:
        print("By complexity score (exec / views_with_charts):")
        for key in sorted(complexity_pairs.keys(), key=_complexity_sort_key):
            pairs = complexity_pairs[key]
            c_ok = sum(1 for s, _ in pairs if s)
            vc_ok = sum(1 for _, v in pairs if v)
            c_n = len(pairs)
            print(f"  complexity {key!r}: exec {c_ok}/{c_n}, charts {vc_ok}/{c_n}")

    return 0 if failure_count == 0 else 1


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--csv-eval-worker":
        sys.exit(_csv_eval_worker_main())
    sys.exit(main())
