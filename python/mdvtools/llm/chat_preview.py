"""
Format subprocess stdout for display in the chat window (truncation, cleanup).
"""

from __future__ import annotations

import re
from typing import Optional

# Defaults aligned with product: large enough for tabular previews, bounded for UI/socket size.
_DEFAULT_MAX_CHARS = 10_000
_DEFAULT_MAX_LINES = 80


def format_stdout_for_chat(
    stdout: Optional[str],
    *,
    max_chars: int = _DEFAULT_MAX_CHARS,
    max_lines: int = _DEFAULT_MAX_LINES,
) -> Optional[str]:
    """
    Normalize stdout into markdown-friendly text for chat, or None if empty.

    Collapses long runs of blank lines, then applies line and character caps with a footer note.
    """
    if stdout is None:
        return None
    text = stdout.replace("\r\n", "\n").strip()
    if not text:
        return None

    text = re.sub(r"\n{3,}", "\n\n", text)

    raw_lines = text.splitlines()
    omitted_lines = 0
    if len(raw_lines) > max_lines:
        omitted_lines = len(raw_lines) - max_lines
        text = "\n".join(raw_lines[:max_lines])
    else:
        text = "\n".join(raw_lines)

    omitted_chars = 0
    if len(text) > max_chars:
        omitted_chars = len(text) - max_chars
        text = text[:max_chars].rstrip()

    if omitted_lines or omitted_chars:
        parts = []
        if omitted_lines:
            parts.append(f"{omitted_lines} lines omitted")
        if omitted_chars:
            parts.append(f"{omitted_chars} characters omitted")
        text = text + "\n\n_" + "; ".join(parts) + "._"

    return text
