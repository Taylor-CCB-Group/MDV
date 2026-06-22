"""Tests for the optional-dependency guard helper (``mdvtools._optional``).

The ``app`` extra (database/server, chat/LLM, auth) is gated behind
``require_extra`` so a slim ``pip install mdvtools`` fails with an actionable
message instead of a bare ``ModuleNotFoundError``. These tests pin both the
helper's behaviour and the *contract* that every heavy leaf module guards itself
— so a removed guard, or a new heavy module added without one, is caught here
rather than by a confused slim-install user.
"""
import re
from pathlib import Path

import pytest

from mdvtools._optional import require_extra, _available

# The leaf modules that import the `app` stack at top level, with the import
# name each one probes. Kept in sync with mdvtools/_optional.py usage.
GUARDED_MODULES = {
    "llm/langchain_mdv.py": "langchain_openai",
    "llm/chatlog.py": "langchain_core",
    "llm/github_utils.py": "nbformat",
    "dbutils/dbmodels.py": "flask_sqlalchemy",
    "dbutils/mdv_server_app.py": "sqlalchemy",
    "auth/auth0_provider.py": "authlib",
}

_MISSING = "definitely_not_installed_pkg_xyz"


def test_available_true_for_stdlib():
    assert _available("json") is True


def test_available_false_for_missing():
    assert _available(_MISSING) is False


def test_require_extra_passes_when_present():
    # Should not raise: stdlib is always importable.
    require_extra("app", "json")


def test_require_extra_raises_with_actionable_message():
    with pytest.raises(ModuleNotFoundError) as exc:
        require_extra("app", _MISSING)
    msg = str(exc.value)
    assert 'pip install "mdvtools[app]"' in msg
    assert _MISSING in msg  # names the missing package


def test_require_extra_reports_all_missing():
    with pytest.raises(ModuleNotFoundError) as exc:
        require_extra("app", "json", _MISSING, "os")
    msg = str(exc.value)
    assert _MISSING in msg
    assert "json" not in msg  # present packages are not listed as missing


@pytest.mark.parametrize("rel_path, probe", sorted(GUARDED_MODULES.items()))
def test_heavy_leaf_modules_are_guarded(rel_path, probe):
    """Each heavy leaf must call require_extra('app', <probe>) before its
    optional imports. Static source check so it holds even in a full install."""
    src = (Path(__file__).resolve().parents[1] / rel_path).read_text()
    assert f'require_extra("app", "{probe}")' in src, (
        f"{rel_path} is missing its app-extra guard for {probe!r}"
    )
    guard_pos = src.index("require_extra(")
    # The guard must precede the real top-level import it protects (a line
    # starting `import <probe>` or `from <probe>` — not the guard call itself).
    import_match = re.search(rf"^(?:import|from) {re.escape(probe)}\b", src, re.MULTILINE)
    assert import_match is not None, f"{rel_path}: no top-level import of {probe!r} found"
    assert guard_pos < import_match.start(), f"{rel_path}: guard must run before importing {probe!r}"
