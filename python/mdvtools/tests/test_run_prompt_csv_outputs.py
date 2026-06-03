"""Tests for slim results CSV + details JSONL split in run_prompt_csv."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

_PY_ROOT = Path(__file__).resolve().parents[2]
if str(_PY_ROOT) not in sys.path:
    sys.path.insert(0, str(_PY_ROOT))

from mdvtools.test_projects.chatmdv_csv_eval import run_prompt_csv


def _write_input_csv(path: Path) -> None:
    path.write_text(
        "Path,Dataset,Question for evaluation,Complexity score\n"
        "/tmp/proj,DS,Test prompt one,1\n",
        encoding="utf-8",
    )


def test_run_prompt_csv_slim_results_and_details_jsonl(tmp_path, monkeypatch):
    input_csv = tmp_path / "prompts.csv"
    out_dir = tmp_path / "report"
    _write_input_csv(input_csv)

    fake_verification = "## What you can verify\n- chart ok\n"

    def fake_run_chat_once(*, project_path, prompt, output_dir=None, view_name=None):
        del project_path, prompt, output_dir, view_name
        return {
            "success": True,
            "message": "Success",
            "captured_output": "Block 'b14' took 1.0 seconds\n",
            "duration_seconds": 3.5,
            "view_name": "test-view",
            "view_snapshot_present": True,
            "chart_count": 1,
            "verification": fake_verification,
            "needs_refresh": False,
            "debug_output_dir": "",
        }, 0

    monkeypatch.setattr(
        "mdvtools.llm.chat_cli.run_chat_once",
        fake_run_chat_once,
    )

    argv = [
        "run_prompt_csv",
        "--csv",
        str(input_csv),
        "--output-dir",
        str(out_dir),
        "--in-process-rows",
    ]
    monkeypatch.setattr(sys, "argv", argv)
    exit_code = run_prompt_csv.main()
    assert exit_code == 0

    results_files = list(out_dir.glob("results_*.csv"))
    assert len(results_files) == 1
    results_path = results_files[0]

    with results_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    assert len(rows) == 1
    excluded = {
        "failure_reason",
        "captured_output_excerpt",
        "verification",
        "run_benchmark_rows_total",
        "run_benchmark_summary_json_mirror",
    }
    for col in excluded:
        assert col not in fieldnames

    assert "exec_success" in fieldnames
    assert "details_jsonl_path" in fieldnames
    assert rows[0]["exec_success"] == "True"
    assert rows[0]["view_name"] == "test-view"

    details_files = list(out_dir.glob("details_*.jsonl"))
    assert len(details_files) == 1
    detail_lines = details_files[0].read_text(encoding="utf-8").strip().splitlines()
    assert len(detail_lines) == 1
    detail = json.loads(detail_lines[0])
    assert detail["verification"] == fake_verification
    assert detail["row_index"] == 1
    assert detail["exec_success"] is True

    summary_files = list(out_dir.glob("benchmark_summary_*.json"))
    assert len(summary_files) == 1
    summary = json.loads(summary_files[0].read_text(encoding="utf-8"))
    assert summary["run_id"]
    assert summary["details_jsonl"] == str(details_files[0])
    assert summary["totals"]["rows"] == 1
    assert "run_benchmark_summary_json_mirror" not in summary


def test_run_prompt_csv_failure_details_in_jsonl_not_csv(tmp_path, monkeypatch):
    input_csv = tmp_path / "prompts.csv"
    out_dir = tmp_path / "report"
    _write_input_csv(input_csv)

    def fake_run_chat_once(*, project_path, prompt, output_dir=None, view_name=None):
        del project_path, prompt, output_dir, view_name
        return {
            "success": False,
            "message": "ERROR: Code execution failed: traceback here",
            "captured_output": "stderr tail",
            "duration_seconds": 1.0,
            "view_name": "",
            "view_snapshot_present": False,
            "chart_count": 0,
            "verification": "",
            "needs_refresh": False,
            "debug_output_dir": "",
        }, 1

    monkeypatch.setattr(
        "mdvtools.llm.chat_cli.run_chat_once",
        fake_run_chat_once,
    )

    argv = [
        "run_prompt_csv",
        "--csv",
        str(input_csv),
        "--output-dir",
        str(out_dir),
        "--in-process-rows",
    ]
    monkeypatch.setattr(sys, "argv", argv)
    exit_code = run_prompt_csv.main()
    assert exit_code == 1

    results_path = next(out_dir.glob("results_*.csv"))
    with results_path.open(newline="", encoding="utf-8") as f:
        row = next(csv.DictReader(f))
    assert "failure_reason" not in row
    assert row["exec_success"] == "False"

    detail = json.loads(next(out_dir.glob("details_*.jsonl")).read_text(encoding="utf-8").strip())
    assert "traceback" in detail["failure_reason"]
    assert detail["captured_output_excerpt"]

    failures_path = next(out_dir.glob("failures_*.jsonl"))
    failure = json.loads(failures_path.read_text(encoding="utf-8").strip())
    assert failure["failure_reason"] == detail["failure_reason"]
