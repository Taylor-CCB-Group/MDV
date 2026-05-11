"""Tests for chat stdout preview formatting."""

from mdvtools.llm.chat_preview import format_stdout_for_chat


def test_empty_and_none():
    assert format_stdout_for_chat(None) is None
    assert format_stdout_for_chat("") is None
    assert format_stdout_for_chat("   \n  ") is None


def test_small_output_unchanged():
    s = "a\nb\nc"
    assert format_stdout_for_chat(s) == s


def test_collapses_excessive_blank_lines():
    s = "x\n\n\n\ny"
    out = format_stdout_for_chat(s)
    assert out is not None
    assert "\n\n\n" not in out


def test_truncates_lines():
    lines = [f"line{i}" for i in range(100)]
    text = "\n".join(lines)
    out = format_stdout_for_chat(text, max_lines=10, max_chars=100_000)
    assert out is not None
    assert "lines omitted" in out
    assert "line0" in out
    assert "line99" not in out


def test_truncates_chars_after_lines():
    text = "x" * 5000
    out = format_stdout_for_chat(text, max_lines=500, max_chars=100)
    assert out is not None
    assert "characters omitted" in out
    assert len(out) <= 200
