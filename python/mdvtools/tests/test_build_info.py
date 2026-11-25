import pytest

from mdvtools import build_info

COMMIT = "dae956a1a976df15761b7360742f9d390ad42dde"

# PEP 610 record written by pip/uv for `mdvtools @ git+https://...#subdirectory=python`
GIT_INSTALL = {
    "url": "https://github.com/Taylor-CCB-Group/MDV",
    "vcs_info": {"vcs": "git", "commit_id": COMMIT, "requested_revision": COMMIT},
    "subdirectory": "python",
}


@pytest.fixture
def no_env_or_checkout(monkeypatch):
    """Simulate an installed copy: no build env vars and no enclosing checkout."""
    monkeypatch.delenv("GIT_COMMIT_HASH", raising=False)
    monkeypatch.setattr(build_info, "_get_from_git", lambda: None)


def test_git_install_reports_its_commit(monkeypatch, no_env_or_checkout):
    monkeypatch.setattr(build_info, "read_direct_url", lambda: GIT_INSTALL)

    info = build_info.get_build_info()

    assert info["git_commit_hash"] == COMMIT
    assert info["repository"] == "https://github.com/Taylor-CCB-Group/MDV (python)"
    assert info["source"] == "git install"
    markdown = build_info.build_info_to_markdown(info)
    assert f"`{COMMIT}`" in markdown
    assert "not available" not in markdown


@pytest.mark.parametrize(
    "record",
    [
        None,
        {"url": "file:///src/MDV/python", "dir_info": {"editable": True}},
        {"url": "https://example.org/mdvtools.tar.gz", "archive_info": {}},
    ],
)
def test_installs_without_a_commit_stay_unknown(monkeypatch, no_env_or_checkout, record):
    monkeypatch.setattr(build_info, "read_direct_url", lambda: record)

    assert build_info.get_build_info()["source"] == "unknown"


def test_checkout_wins_over_direct_url(monkeypatch):
    monkeypatch.delenv("GIT_COMMIT_HASH", raising=False)
    checkout = {"git_commit_hash": "abc", "source": "git"}
    monkeypatch.setattr(build_info, "_get_from_git", lambda: checkout)
    monkeypatch.setattr(build_info, "read_direct_url", lambda: GIT_INSTALL)

    assert build_info.get_build_info() is checkout
