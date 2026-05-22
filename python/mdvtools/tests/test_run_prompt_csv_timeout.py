"""Per-row subprocess timeout in run_prompt_csv."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

_PY_ROOT = Path(__file__).resolve().parents[2]
if str(_PY_ROOT) not in sys.path:
    sys.path.insert(0, str(_PY_ROOT))

from mdvtools.test_projects.chatmdv_csv_eval import run_prompt_csv


def test_run_chat_row_subprocess_timeout_returns_failed_result():
    exc = subprocess.TimeoutExpired(
        cmd=["python"],
        timeout=10.0,
        output=b"partial stdout",
        stderr=b"partial stderr",
    )
    with patch.object(run_prompt_csv.subprocess, "run", side_effect=exc):
        result, exit_code, term_tag = run_prompt_csv._run_chat_row_subprocess(
            project_path="/tmp/proj",
            prompt="test",
            output_dir=None,
            timeout_seconds=10.0,
        )
    assert term_tag == "timeout"
    assert exit_code == 5
    assert result["success"] is False
    assert "timed out after 10" in str(result["message"])
    assert "partial stdout" in str(result["captured_output"])
