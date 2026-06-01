"""Helpers for optional dependency groups (PEP 621 extras).

The database/server, chat/LLM and auth stack ships behind the single ``app``
extra (``pip install "mdvtools[app]"``). A bare ``pip install mdvtools`` is the
slim core (project manipulation, viewing, ``MDVProject.serve()``).

Modules that import the ``app`` stack at top level call :func:`require_extra`
*before* those imports, so a slim install fails with an actionable message
instead of a bare ``ModuleNotFoundError`` from a transitive package name the
user has never heard of.
"""
from __future__ import annotations

import importlib.util


def _available(package: str) -> bool:
    """Return True if ``package`` can be imported, without importing it."""
    try:
        return importlib.util.find_spec(package) is not None
    except (ImportError, ValueError):
        # A parent package is itself missing/broken — treat as unavailable.
        return False


def require_extra(extra: str, *packages: str) -> None:
    """Raise a friendly error if any of ``packages`` is not importable.

    Call at the top of a module, above its optional top-level imports::

        from mdvtools._optional import require_extra
        require_extra("app", "langchain_openai")

        from langchain_openai import ChatOpenAI  # safe below the guard

    ``packages`` are *import* names to probe (e.g. ``"flask_sqlalchemy"``),
    chosen to be present only when ``extra`` is installed.
    """
    missing = [p for p in packages if not _available(p)]
    if missing:
        raise ModuleNotFoundError(
            f"This feature requires the optional '{extra}' dependencies "
            f"(missing: {', '.join(missing)}).\n"
            f'Install them with:  pip install "mdvtools[{extra}]"'
        )
