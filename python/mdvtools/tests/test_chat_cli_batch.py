import csv
from pathlib import Path

from mdvtools.llm import chat_cli
from mdvtools.mdvproject import MDVProject


def _make_project(tmp_path):
    project_dir = tmp_path / "proj"
    project = MDVProject(str(project_dir))
    project.datasources = [{"name": "cells", "columns": [{"field": "cell_id"}]}]
    return project


def test_validate_mode_args_rules():
    parser = chat_cli._build_parser()
    args = parser.parse_args(["--project", "/tmp/x", "--prompt", "a"])
    chat_cli._validate_mode_args(args)

    args = parser.parse_args(["--project", "/tmp/x", "--prompt-file", "/tmp/p.txt"])
    try:
        chat_cli._validate_mode_args(args)
        assert False, "expected ValueError"
    except ValueError as exc:
        assert "--csv-log is required" in str(exc)


def test_batch_prompts_continue_on_error_and_write_csv(tmp_path, monkeypatch):
    project = _make_project(tmp_path)
    prompt_file = tmp_path / "prompts.txt"
    prompt_file.write_text("first prompt\nsecond prompt\n", encoding="utf-8")
    csv_log = tmp_path / "bench.csv"
    output_dir = tmp_path / "runs"

    calls = {"n": 0}

    def fake_run_chat_once(*, project_path, prompt, output_dir=None, view_name=None):
        del project_path, view_name
        calls["n"] += 1
        run_out = Path(output_dir)
        run_out.mkdir(parents=True, exist_ok=True)
        (run_out / "generated_code.py").write_text("print('ok')\n", encoding="utf-8")
        (run_out / "result.json").write_text("{}", encoding="utf-8")
        (run_out / "view_snapshot.json").write_text("{}", encoding="utf-8")
        if prompt.startswith("second"):
            return {
                "success": False,
                "project_path": str(project.dir),
                "views_file": str(Path(project.dir) / "views.json"),
                "view_name": None,
                "message": "boom",
                "debug_output_dir": str(run_out.resolve()),
                "block_timings": {"block_b14_execute_code_s": 1.2},
                "captured_output": "",
                "duration_seconds": 2.0,
            }, 2
        return {
            "success": True,
            "project_path": str(project.dir),
            "views_file": str(Path(project.dir) / "views.json"),
            "view_name": "v1",
            "message": "Success",
            "debug_output_dir": str(run_out.resolve()),
            "block_timings": {"block_b12_rag_chain_s": 3.4},
            "captured_output": "",
            "duration_seconds": 4.0,
        }, 0

    monkeypatch.setattr(chat_cli, "run_chat_once", fake_run_chat_once)
    rows, exit_code = chat_cli.run_batch_prompts(
        project_path=str(project.dir),
        prompt_file=str(prompt_file),
        csv_log=str(csv_log),
        base_url="http://localhost:5050",
        output_dir=str(output_dir),
    )

    assert calls["n"] == 2
    assert exit_code == 1
    assert len(rows) == 2
    assert rows[0]["status"] == "success"
    assert rows[1]["status"] == "failure"
    assert rows[0]["view_url"] == "http://localhost:5050/project/proj/view/v1"
    assert rows[1]["view_url"] == ""

    with csv_log.open("r", encoding="utf-8") as f:
        records = list(csv.DictReader(f))
    assert len(records) == 2
    assert "block_b12_rag_chain_s" in records[0]
    assert "block_b14_execute_code_s" in records[0]
    assert records[1]["error_message"] == "boom"


def test_parse_block_timings():
    text = "Block 'b12: RAG chain' took 2.5000 seconds\nBlock 'b14: Execute code' took 0.1250 seconds"
    timings = chat_cli._parse_block_timings(text)
    assert timings["block_b12_rag_chain_s"] == 2.5
    assert timings["block_b14_execute_code_s"] == 0.125
