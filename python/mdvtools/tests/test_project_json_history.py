"""Tests for the standalone project JSON history CLI."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "project_json_history.py"
_spec = importlib.util.spec_from_file_location("project_json_history", _SCRIPT)
if _spec is None or _spec.loader is None:
    raise RuntimeError(f"could not load {_SCRIPT}")
_history = importlib.util.module_from_spec(_spec)
sys.modules["project_json_history"] = _history
_spec.loader.exec_module(_history)

discover_projects = getattr(_history, "discover_projects")
main = getattr(_history, "main")


def _write_project(
    path: Path,
    datasources: bytes,
    views: bytes | None = None,
    state: bytes | None = None,
) -> None:
    path.mkdir(parents=True, exist_ok=True)
    (path / "datasources.json").write_bytes(datasources)
    if views is not None:
        (path / "views.json").write_bytes(views)
    if state is not None:
        (path / "state.json").write_bytes(state)


def _parse_snapshot(text: str) -> tuple[str, dict[str, str]]:
    lines = [line for line in text.splitlines() if line]
    assert lines[0].startswith("run ")
    run = lines[0].removeprefix("run ")
    mapping: dict[str, str] = {}
    for line in lines[1:]:
        project, rev = line.rsplit(" ", 1)
        mapping[project] = rev
    return run, mapping


def _manifest(project: Path, rev: str) -> dict[str, object]:
    payload = json.loads((project / "history" / rev / "manifest.json").read_text(encoding="utf-8"))
    assert isinstance(payload, dict)
    return payload


def _head(project: Path) -> str:
    return (project / "history" / "HEAD").read_text(encoding="utf-8").strip()


def test_snapshot_versions_each_project_once_and_shares_a_run(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    liver = tmp_path / "liver"
    pbmc = tmp_path / "pbmc3k"
    views = b'{"default": {"viewImage": "data:image/png;base64,AAAA"}}\n'
    _write_project(liver, b'[{"name": "cells"}]\n', views, b'{"all_views": ["default"]}\n')
    _write_project(pbmc, b'[{"name": "cells", "columns": [{"field": "Cell Type"}]}]\n', views, b"{}\n")
    nested = pbmc / "nested"
    nested.mkdir()
    (nested / "datasources.json").write_bytes(b'[{"name": "ignored"}]\n')
    (pbmc / "datafile.h5").write_bytes(b"hdf5-bytes")

    assert main(["snapshot", str(tmp_path), "--message", "before column rename"]) == 0
    run, mapping = _parse_snapshot(capsys.readouterr().out)

    assert set(mapping) == {str(liver.resolve()), str(pbmc.resolve())}
    assert mapping[str(liver.resolve())] != mapping[str(pbmc.resolve())]
    projects, errors = discover_projects([tmp_path])
    assert errors == []
    assert projects == [liver.resolve(), pbmc.resolve()]

    for project in (liver, pbmc):
        rev = mapping[str(project.resolve())]
        manifest = _manifest(project, rev)
        assert manifest["run"] == run
        assert manifest["parent"] is None
        assert manifest["message"] == "before column rename"
        assert manifest["files"] == ["datasources.json", "views.json", "state.json"]
        assert _head(project) == rev
        assert (project / "history" / rev / "views.json").read_bytes() == views
        assert not (project / "nested" / "history").exists()
        assert not (project / "history" / rev / "history").exists()
    assert (pbmc / "datafile.h5").read_bytes() == b"hdf5-bytes"


def test_second_snapshot_records_parent_and_a_new_run(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    project = tmp_path / "pbmc3k"
    _write_project(project, b'{"v": 1}\n', b'{"view": 1}\n')
    assert main(["snapshot", str(project), "--message", "before column rename"]) == 0
    first_run, first = _parse_snapshot(capsys.readouterr().out)
    first_id = first[str(project.resolve())]

    (project / "datasources.json").write_bytes(b'{"v": 2}\n')
    assert main(["snapshot", str(project), "--message", "after column rename"]) == 0
    second_run, second = _parse_snapshot(capsys.readouterr().out)
    second_id = second[str(project.resolve())]

    assert first_run != second_run
    assert _head(project) == second_id
    manifest = _manifest(project, second_id)
    assert manifest["parent"] == first_id
    assert manifest["run"] == second_run
    assert (project / "history" / second_id / "datasources.json").read_bytes() == b'{"v": 2}\n'
    assert (project / "history" / first_id / "datasources.json").read_bytes() == b'{"v": 1}\n'


def test_project_with_only_datasources_json(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    project = tmp_path / "notes-only"
    _write_project(project, b'[{"name": "cells"}]\n')
    assert main(["snapshot", str(project)]) == 0
    _run, mapping = _parse_snapshot(capsys.readouterr().out)
    rev = mapping[str(project.resolve())]
    assert _manifest(project, rev)["files"] == ["datasources.json"]
    assert not (project / "history" / rev / "views.json").exists()


def test_log_one_project_and_a_tree(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    liver = tmp_path / "liver"
    pbmc = tmp_path / "pbmc3k"
    _write_project(liver, b"{}\n")
    _write_project(pbmc, b"{}\n")
    assert main(["snapshot", str(tmp_path), "--message", "before column rename"]) == 0
    capsys.readouterr()
    (pbmc / "datasources.json").write_bytes(b'{"changed": true}\n')
    assert main(["snapshot", str(pbmc), "--message", "after column rename"]) == 0
    capsys.readouterr()

    assert main(["log", str(pbmc)]) == 0
    one = capsys.readouterr().out
    assert str(pbmc.resolve()) not in one
    assert "parent=-" in one
    assert "run=" in one
    assert "before column rename" in one
    assert "after column rename" in one
    assert one.index("after column rename") < one.index("before column rename")

    assert main(["log", str(tmp_path)]) == 0
    tree = capsys.readouterr().out
    assert str(liver.resolve()) in tree
    assert str(pbmc.resolve()) in tree


def test_restore_one_revision_puts_bytes_back_and_appends_history(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    project = tmp_path / "pbmc3k"
    original_ds = b'[{"field": "Cell Type"}]\n'
    original_views = b'{"param": "Cell Type", "viewImage": "data:image/png;base64,AAAA"}\n'
    _write_project(project, original_ds, original_views, b'{"all_views": ["default"]}\n')
    (project / "datafile.h5").write_bytes(b"leave-me")
    assert main(["snapshot", str(project), "--message", "before column rename"]) == 0
    _run, first = _parse_snapshot(capsys.readouterr().out)
    first_id = first[str(project.resolve())]

    (project / "datasources.json").write_bytes(b'[{"field": "cell_type"}]\n')
    (project / "views.json").write_bytes(b'{"param": "cell_type"}\n')
    assert main(["snapshot", str(project), "--message", "after column rename"]) == 0
    capsys.readouterr()
    pre_restore_head = _head(project)

    assert main(["restore", str(project), first_id]) == 0
    restored = capsys.readouterr().out.strip()
    assert restored == f"{project.resolve()} {first_id}"
    assert (project / "datasources.json").read_bytes() == original_ds
    assert (project / "views.json").read_bytes() == original_views
    assert (project / "datafile.h5").read_bytes() == b"leave-me"
    new_head = _head(project)
    assert new_head != pre_restore_head
    manifest = _manifest(project, new_head)
    assert manifest["parent"] == pre_restore_head
    assert manifest["message"] == f"restore {first_id}"
    assert (project / "history" / new_head / "views.json").read_bytes() == original_views


def test_restore_run_restores_matching_projects_and_skips_the_rest(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
):
    liver = tmp_path / "liver"
    pbmc = tmp_path / "pbmc3k"
    _write_project(liver, b'{"name": "liver"}\n', b'{"v": 1}\n')
    _write_project(pbmc, b'{"name": "pbmc"}\n', b'{"v": 1}\n')
    assert main(["snapshot", str(tmp_path), "--message", "before column rename"]) == 0
    run, mapping = _parse_snapshot(capsys.readouterr().out)
    (liver / "datasources.json").write_bytes(b'{"name": "changed-liver"}\n')
    (pbmc / "datasources.json").write_bytes(b'{"name": "changed-pbmc"}\n')

    added = tmp_path / "added-later"
    _write_project(added, b'{"name": "added"}\n')

    assert main(["restore", str(tmp_path), "--run", run]) == 0
    captured = capsys.readouterr()
    assert f"{liver.resolve()} {mapping[str(liver.resolve())]}" in captured.out
    assert f"{pbmc.resolve()} {mapping[str(pbmc.resolve())]}" in captured.out
    assert f"skip {added.resolve()}: no revision for run {run}" in captured.err
    assert (liver / "datasources.json").read_bytes() == b'{"name": "liver"}\n'
    assert (pbmc / "datasources.json").read_bytes() == b'{"name": "pbmc"}\n'
    assert (added / "datasources.json").read_bytes() == b'{"name": "added"}\n'
    assert not (added / "history").exists()
    liver_head = _manifest(liver, _head(liver))
    pbmc_head = _manifest(pbmc, _head(pbmc))
    assert liver_head["run"] == pbmc_head["run"]
    assert liver_head["run"] != run
    assert liver_head["message"] == f"restore {mapping[str(liver.resolve())]}"


def test_restore_failures_do_not_write(tmp_path: Path, capsys: pytest.CaptureFixture[str]):
    liver = tmp_path / "liver"
    pbmc = tmp_path / "pbmc3k"
    _write_project(liver, b'{"name": "liver"}\n')
    _write_project(pbmc, b'{"name": "pbmc"}\n', b'{"view": 1}\n')
    assert main(["snapshot", str(tmp_path), "--message", "before column rename"]) == 0
    run, mapping = _parse_snapshot(capsys.readouterr().out)
    pbmc_rev = mapping[str(pbmc.resolve())]
    (pbmc / "datasources.json").write_bytes(b'{"name": "edited"}\n')
    (pbmc / "views.json").write_bytes(b'{"view": 2}\n')
    edited_ds = (pbmc / "datasources.json").read_bytes()
    edited_views = (pbmc / "views.json").read_bytes()
    head_before = _head(pbmc)

    assert main(["restore", str(tmp_path), pbmc_rev]) == 1
    err = capsys.readouterr().err
    assert "--run" in err
    assert (pbmc / "datasources.json").read_bytes() == edited_ds
    assert _head(pbmc) == head_before

    assert main(["restore", str(pbmc), "20260927T080000Z-ffffff"]) == 1
    capsys.readouterr()
    assert (pbmc / "datasources.json").read_bytes() == edited_ds
    assert _head(pbmc) == head_before

    assert main(["restore", str(tmp_path), "--run", "0000000000000000"]) == 1
    missing = capsys.readouterr().err
    assert "no project has run 0000000000000000" in missing
    assert (pbmc / "datasources.json").read_bytes() == edited_ds
    assert _head(pbmc) == head_before

    (pbmc / "history" / pbmc_rev / "datasources.json").unlink()
    assert main(["restore", str(tmp_path), "--run", run]) == 1
    capsys.readouterr()
    assert (pbmc / "datasources.json").read_bytes() == edited_ds
    assert (pbmc / "views.json").read_bytes() == edited_views
    assert (liver / "datasources.json").read_bytes() == b'{"name": "liver"}\n'
    assert _head(pbmc) == head_before
    assert _head(liver) == mapping[str(liver.resolve())]
