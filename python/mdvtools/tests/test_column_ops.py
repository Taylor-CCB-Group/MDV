import pandas as pd
import pytest
from click.testing import CliRunner

from mdvtools.cli import cli
from mdvtools.mdvproject import MDVProject


def _build_project(tmp_path):
    project = MDVProject(str(tmp_path), delete_existing=True)
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    project.add_datasource("cells", df)
    return project


def test_drop_columns_soft_keeps_h5_dataset(tmp_path):
    project = _build_project(tmp_path)
    report = project.drop_columns("cells", ["a"], hard=False, strict=False)
    assert report["dropped"] == ["a"]
    metadata = project.get_datasource_metadata("cells")
    fields = {col["field"] for col in metadata["columns"]}
    assert "a" not in fields
    h5 = project._get_h5_handle(read_only=True)
    try:
        assert "a" in h5["cells"]
    finally:
        h5.close()


def test_drop_columns_hard_removes_h5_dataset(tmp_path):
    project = _build_project(tmp_path)
    report = project.drop_columns("cells", ["a"], hard=True, strict=False)
    assert report["hard"] is True
    h5 = project._get_h5_handle(read_only=True)
    try:
        assert "a" not in h5["cells"]
    finally:
        h5.close()


def test_rename_columns_display_name_only(tmp_path):
    project = _build_project(tmp_path)
    report = project.rename_columns("cells", {"a": "Gene A"}, rename_field=False)
    assert report["renamed"] == {"a": "Gene A"}
    metadata = project.get_datasource_metadata("cells")
    target = next(col for col in metadata["columns"] if col["field"] == "a")
    assert target["name"] == "Gene A"
    h5 = project._get_h5_handle(read_only=True)
    try:
        assert "a" in h5["cells"]
    finally:
        h5.close()


def test_rename_field_requires_hard(tmp_path):
    project = _build_project(tmp_path)
    with pytest.raises(ValueError, match="requires hard=True"):
        project.rename_columns("cells", {"a": "z"}, rename_field=True, hard=False)


def test_hard_field_rename_updates_h5_and_references(tmp_path):
    project = _build_project(tmp_path)
    project.views = {"main": {"initialCharts": {"cells": [{"param": "a"}]}}}
    state = project.state
    state["filters"] = {"cells": {"column": "a"}}
    project.state = state
    ds = project.get_datasource_metadata("cells")
    ds["interactions"] = {"pivot_column": "a"}
    project.set_datasource_metadata(ds)

    project.rename_columns("cells", {"a": "z"}, rename_field=True, hard=True)

    metadata = project.get_datasource_metadata("cells")
    fields = {col["field"] for col in metadata["columns"]}
    assert "z" in fields
    assert "a" not in fields
    assert metadata["interactions"]["pivot_column"] == "z"
    assert project.views["main"]["initialCharts"]["cells"][0]["param"] == "z"
    assert project.state["filters"]["cells"]["column"] == "z"
    h5 = project._get_h5_handle(read_only=True)
    try:
        assert "z" in h5["cells"]
        assert "a" not in h5["cells"]
    finally:
        h5.close()


def test_cli_field_id_requires_hard(tmp_path):
    project = _build_project(tmp_path)
    assert project.get_datasource_metadata("cells")
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "rename-columns",
            str(tmp_path),
            "--datasource",
            "cells",
            "--map",
            "a:z",
            "--field-id",
        ],
    )
    assert result.exit_code != 0
    assert "--field-id requires --hard" in result.output


def test_cli_restore_backup_does_not_require_datasource(tmp_path):
    project = _build_project(tmp_path)
    # create a backup
    project.drop_columns("cells", ["a"], strict=False, cleanup=True, backup=True)
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "drop-columns",
            str(tmp_path),
            "--restore-backup",
        ],
    )
    assert result.exit_code == 0
    assert "restoredFiles" in result.output


def test_drop_columns_cleanup_removes_from_views_and_creates_backup(tmp_path):
    project = _build_project(tmp_path)
    project.views = {
        "Overview Screen": {
            "initialCharts": {
                "cells": [
                    {
                        "type": "table_chart_react",
                        "id": "NpEPWW",
                        "param": ["a", "b"],
                    }
                ]
            }
        }
    }
    report = project.drop_columns(
        "cells",
        ["a"],
        strict=True,
        hard=False,
        cleanup=True,
        backup=True,
    )
    assert report["cleanup"] is True
    assert report["backup"]["backupDir"] is not None
    # column 'a' removed from table chart params
    assert project.views["Overview Screen"]["initialCharts"]["cells"][0]["param"] == ["b"]
    # backup directory created and contains views.json
    backup_dir = report["backup"]["backupDir"]
    assert backup_dir
    assert "views.json" in report["backup"]["backedUp"]


def test_drop_columns_cleanup_can_fix_already_missing_field(tmp_path):
    project = _build_project(tmp_path)
    # simulate already-missing column metadata (like DOK3) but still referenced in views
    project.drop_columns("cells", ["a"], strict=False, hard=False)
    project.views = {
        "Overview Screen": {
            "initialCharts": {
                "cells": [
                    {
                        "type": "table_chart_react",
                        "id": "NpEPWW",
                        "param": ["a", "b"],
                    }
                ]
            }
        }
    }
    # strict=True should still allow cleanup to remove the missing field from views
    report = project.drop_columns(
        "cells",
        ["a"],
        strict=True,
        hard=False,
        cleanup=True,
        backup=False,
    )
    assert report["missing"] == ["a"]
    assert project.views["Overview Screen"]["initialCharts"]["cells"][0]["param"] == ["b"]


def test_restore_json_backup_latest(tmp_path):
    project = _build_project(tmp_path)
    project.views = {"main": {"initialCharts": {"cells": [{"type": "table_chart_react", "param": ["a"]}]}}}
    project.drop_columns("cells", ["a"], strict=False, cleanup=True, backup=True)
    # mutate views after backup
    project.views = {"main": {"initialCharts": {"cells": [{"type": "table_chart_react", "param": ["b"]}]}}}
    report = project.restore_json_backup()
    assert "views.json" in report["restoredFiles"]
    assert project.views["main"]["initialCharts"]["cells"][0]["param"] == ["a"]


def test_restore_json_backup_specific_timestamp(tmp_path):
    project = _build_project(tmp_path)
    project.views = {"v1": {"initialCharts": {"cells": [{"type": "table_chart_react", "param": ["a"]}]}}}
    r1 = project.drop_columns("cells", ["a"], strict=False, cleanup=True, backup=True)
    ts1 = r1["backup"]["backupDir"].split("/")[-1]
    project.views = {"v2": {"initialCharts": {"cells": [{"type": "table_chart_react", "param": ["b"]}]}}}
    # ensure second backup has a different timestamp directory
    import time
    time.sleep(1.05)
    project.drop_columns("cells", ["b"], strict=False, cleanup=True, backup=True)
    # now restore first backup explicitly
    project.restore_json_backup(timestamp=ts1)
    assert "v1" in project.views


def test_restore_json_backup_missing_raises(tmp_path):
    project = _build_project(tmp_path)
    with pytest.raises(FileNotFoundError):
        project.restore_json_backup()
