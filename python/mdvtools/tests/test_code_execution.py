"""Tests for ChatMDV subprocess wrapper (mdvtools.llm.code_execution)."""

from unittest.mock import patch

import subprocess as sp

from mdvtools.llm import code_execution


def test_run_subprocess_calledprocesserror_includes_diagnostic():
    e = sp.CalledProcessError(1, ["python", "x.py"], output=b"hello", stderr=b"errline")
    with patch.object(code_execution.subprocess, "run", side_effect=e):
        out, err = code_execution.run_subprocess(["python", "x.py"])
    assert out is None
    assert "returncode=1" in err
    assert "errline" in err


def test_run_subprocess_sigkill_style_empty_streams():
    e = sp.CalledProcessError(-9, ["python", "x.py"], output=b"", stderr=b"")
    with patch.object(code_execution.subprocess, "run", side_effect=e):
        out, err = code_execution.run_subprocess(["python", "x.py"])
    assert out is None
    assert "returncode=-9" in err


def test_execute_code_failure_includes_source_for_syntax_error():
    bad = "print('ok'\n"  # SyntaxError: unmatched '('
    ok, stdout, stderr = code_execution.execute_code(bad, log=lambda *_a, **_k: None)
    assert ok is False
    assert stdout is None
    assert "Failed Python source" in stderr
    assert bad in stderr
    assert "SyntaxError" in stderr or "unmatched" in stderr.lower()
