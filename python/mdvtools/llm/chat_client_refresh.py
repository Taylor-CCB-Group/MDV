"""
Detect when chat-executed code may have changed project datasources; the browser
should reload so ChartManager receives updated datasources.json metadata.
"""
from __future__ import annotations

import re


def client_needs_refresh_after_chat(code: str | None) -> bool:
    if not code:
        return False
    return bool(
        re.search(r"\badd_datasource\s*\(", code)
        or re.search(r"\badd_datasource_polars\s*\(", code)
        or re.search(r"\bdelete_datasource\s*\(", code)
    )
