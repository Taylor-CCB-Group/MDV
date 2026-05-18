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
or an explicit ``--results`` / ``--failure-log`` path.

Input CSV columns (flexible names)
----------------------------------

- **Path** (or ``path``, ``project_path``): MDV project directory.
- **Question** (or ``Question for evaluation``, ``prompt``): natural-language prompt.
- Optional: **Dataset**, **Complexity score** (or ``complexity``)—copied through and
  included in failure records.

Outputs
-------

- **Results CSV** (default ``results_<stem>_<utc_timestamp>.csv`` under the report directory,
  **outside the repo**, typically ``$TMPDIR/chatmdv_csv_eval/`` — override with ``--output-dir``
  or env ``CHATMDV_CSV_EVAL_OUTPUT_DIR``): original columns plus run metrics:

  - ``exec_success`` — Chat MDV pipeline completed without ``error`` (same as
    ``chat_cli`` ``success``).
  - ``exit_code``, ``duration_seconds``, ``failure_reason`` (empty on success),
    ``captured_output_excerpt`` (tail of captured stdout/stderr on failure).
  - ``view_name`` — name returned after ``ask_question`` (may be empty).
  - ``view_snapshot_present`` — in-memory MDV ``get_view`` returned a snapshot.
  - ``chart_count`` — charts counted under ``view[\"initialCharts\"]``.
  - ``view_has_charts`` — ``exec_success`` and ``chart_count > 0`` (structural
    “built a plotted view”; text-only successes stay False).
  - ``verification``, ``needs_refresh`` — copied from Chat MDV (see
    ``AskQuestionResult``).

  **Replicated run summary (same on every result row):** ``run_id``,
  ``run_completed_at_utc``, ``run_wall_seconds_total``, ``run_cli_limit``,
  ``run_artifacts_base_dir``, paths to this run’s input/results/summary/failure
  logs, row and failure counts, exec and ``view_has_charts`` rates (0–1 and
  percent), compact JSON columns ``run_benchmark_by_complexity_json``,
  ``run_benchmark_totals_json``, and ``run_benchmark_summary_json_mirror`` (full
  aggregate payload matching the JSON file, minified).
- **`benchmark_summary_<stem>_<timestamp>.json`** (always written, same folder as results):

  .. code-block:: json

    {
      "benchmark_summary": "<path to this file>",
      "by_complexity": {"1": {"exec_success_count": N, "rows": N, "view_has_charts_count": N}},
      "failure_log": "<path or null>",
      "input_csv": "<path>",
      "results_csv": "<path>",
      "totals": {
        "exec_success_count": N,
        "exec_success_rate": 0.0,
        "rows": N,
        "view_has_charts_count": N,
        "view_has_charts_rate": 0.0
      }
    }

- **`failures_<stem>_<timestamp>.jsonl`**: written only if at least one row failed;
  one JSON object per failed row (includes full ``captured_output`` plus the view
  fields above).
- **Stderr**: one block per failure as it occurs.
- **Stdout**: summary counts, then a consolidated **Output directory / Files written**
  listing (results CSV, benchmark summary JSON, failures JSONL or why skipped,
  artifacts dir status), complexity breakdown ``exec … / charts …``.

Use ``--artifacts-dir`` to persist per-row ``generated_code.py`` / ``result.json``
(see ``mdvtools.llm.chat_cli``).

Per-row ``exit_code`` matches ``chat_cli.run_chat_once`` (0 success, 1 logical failure,
2 exception). Process exit: **0** if all rows succeeded, **1** if any row failed,
**2** only for bad paths or empty CSV header.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import tempfile
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

_PY_ROOT = Path(__file__).resolve().parents[3]
if str(_PY_ROOT) not in sys.path:
    sys.path.insert(0, str(_PY_ROOT))


EXCERPT_MAX_LEN = 4000


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


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--csv", "-c", required=True, type=Path, help="Benchmark CSV path.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for results, benchmark summary JSON, and failure JSONL (default: "
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


def main() -> int:
    args = _parse_args()
    from mdvtools.llm.chat_cli import run_chat_once
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

    with csv_path.open(newline="", encoding="utf-8") as f:
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
    complexity_pairs: dict[str, list[tuple[bool, bool]]] = defaultdict(list)
    failure_count = 0
    failure_payloads: list[dict[str, object]] = []

    comp_header = _find_header(raw_fieldnames, "complexity")

    for idx, row in enumerate(rows_in, start=1):
        path = _resolve_row_field(row, "path")
        question = _resolve_row_field(row, "question")
        dataset = _resolve_row_field(row, "dataset") or ""
        complexity_val = ((row.get(comp_header) or "").strip() if comp_header else "") or ""

        if not path or not question:
            msg = f"Row {idx}: missing Path or Question (path={path!r}, question={question!r})"
            print(msg, file=sys.stderr)
            out_row = {
                **{k: row.get(k, "") for k in raw_fieldnames},
                "exec_success": False,
                "exit_code": 2,
                "duration_seconds": 0.0,
                "failure_reason": msg,
                "captured_output_excerpt": "",
                "view_name": "",
                "view_snapshot_present": False,
                "chart_count": 0,
                "view_has_charts": False,
                "verification": "",
                "needs_refresh": False,
            }
            result_rows.append(out_row)
            failure_count += 1
            if complexity_val:
                complexity_pairs[complexity_val].append((False, False))
            failure_payloads.append(
                {
                    "row_index": idx,
                    "path": path or "",
                    "dataset": dataset,
                    "question": question or "",
                    "complexity_score": complexity_val,
                    "exit_code": 2,
                    "failure_reason": msg,
                    "duration_seconds": 0.0,
                    "captured_output": "",
                    "view_name": "",
                    "view_snapshot_present": False,
                    "chart_count": 0,
                    "view_has_charts": False,
                    "verification": "",
                    "needs_refresh": False,
                    "debug_output_dir": "",
                }
            )
            continue

        run_out = None
        if base_artifacts is not None:
            run_out = base_artifacts / f"row_{idx:04d}_{_slug(question)}"

        result, exit_code = run_chat_once(
            project_path=path,
            prompt=question,
            output_dir=str(run_out) if run_out else None,
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

        out_row = {
            **{k: row.get(k, "") for k in raw_fieldnames},
            "exec_success": success,
            "exit_code": exit_code,
            "duration_seconds": duration,
            "failure_reason": failure_reason,
            "captured_output_excerpt": excerpt,
            "view_name": view_name,
            "view_snapshot_present": view_snapshot_present,
            "chart_count": chart_count,
            "view_has_charts": view_has_charts,
            "verification": verification,
            "needs_refresh": needs_refresh,
        }
        result_rows.append(out_row)

        if complexity_val:
            complexity_pairs[complexity_val].append((success, view_has_charts))

        if not success:
            failure_count += 1
            q_short = question[:120] + ("..." if len(question) > 120 else "")
            print(
                f"[FAIL] row={idx} exit={exit_code} path={path} q={q_short!r}\n"
                f"       reason: {failure_reason}",
                file=sys.stderr,
            )
            failure_payloads.append(
                {
                    "row_index": idx,
                    "path": path,
                    "dataset": dataset,
                    "question": question,
                    "complexity_score": complexity_val,
                    "exit_code": exit_code,
                    "failure_reason": failure_reason,
                    "duration_seconds": duration,
                    "captured_output": cap,
                    "view_name": view_name,
                    "view_snapshot_present": view_snapshot_present,
                    "chart_count": chart_count,
                    "view_has_charts": view_has_charts,
                    "verification": verification,
                    "needs_refresh": needs_refresh,
                    "debug_output_dir": str(debug_dir) if debug_dir else "",
                }
            )

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

    summary_payload = {
        "input_csv": str(csv_path),
        "results_csv": str(results_path),
        "failure_log": str(failure_log_path) if failure_count > 0 else None,
        "benchmark_summary": str(benchmark_summary_path),
        "totals": {
            "rows": n,
            "exec_success_count": ok,
            "exec_success_rate": round(ok / n, 6) if n else 0.0,
            "view_has_charts_count": view_charts_ok,
            "view_has_charts_rate": round(view_charts_ok / n, 6) if n else 0.0,
        },
        "by_complexity": by_complexity,
    }

    completion_iso = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S+00:00")
    wall_seconds = round(time.perf_counter() - wall_start, 6)

    run_fields = {
        "run_id": utc_now,
        "run_completed_at_utc": completion_iso,
        "run_wall_seconds_total": wall_seconds,
        "run_cli_limit": args.limit if args.limit is not None else "",
        "run_artifacts_base_dir": str(base_artifacts) if base_artifacts else "",
        "run_input_csv_path": str(csv_path),
        "run_results_csv_path": str(results_path),
        "run_benchmark_summary_json_path": str(benchmark_summary_path),
        "run_failure_jsonl_path": str(failure_log_path) if failure_count > 0 else "",
        "run_benchmark_rows_total": n,
        "run_benchmark_failure_count": failure_count,
        "run_benchmark_exec_success_count": ok,
        "run_benchmark_exec_success_rate": round(ok / n, 6) if n else 0.0,
        "run_benchmark_exec_success_pct": round(ok / n * 100.0, 4) if n else 0.0,
        "run_benchmark_view_has_charts_count": view_charts_ok,
        "run_benchmark_view_has_charts_rate": round(view_charts_ok / n, 6) if n else 0.0,
        "run_benchmark_view_has_charts_pct": round(view_charts_ok / n * 100.0, 4) if n else 0.0,
        "run_benchmark_by_complexity_json": json.dumps(
            by_complexity, separators=(",", ":"), sort_keys=True
        ),
        "run_benchmark_totals_json": json.dumps(
            summary_payload["totals"], separators=(",", ":"), sort_keys=True
        ),
        "run_benchmark_summary_json_mirror": json.dumps(
            summary_payload, separators=(",", ":"), sort_keys=True
        ),
    }
    run_aggregate_fields = list(run_fields.keys())
    for row in result_rows:
        row.update(run_fields)

    extra_fields = [
        "exec_success",
        "exit_code",
        "duration_seconds",
        "failure_reason",
        "captured_output_excerpt",
        "view_name",
        "view_snapshot_present",
        "chart_count",
        "view_has_charts",
        "verification",
        "needs_refresh",
    ]
    out_fieldnames = list(raw_fieldnames)
    for f in extra_fields + run_aggregate_fields:
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

    print("Output directory (persisted artifacts):")
    print(f"  {report_dir}")
    print("Files written:")
    print(f"  results_csv               {results_path}")
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
        "During the run: brief [FAIL] lines on stderr."
        + (
            f' Full payloads (captured_output, view fields): see failures_jsonl above ("{failure_log_path.name}").'
            if failure_count
            else " No failures_jsonl written (nothing failed)."
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
    sys.exit(main())
