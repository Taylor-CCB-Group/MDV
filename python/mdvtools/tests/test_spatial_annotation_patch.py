import os

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from mdvtools.mdvproject import MDVProject
from mdvtools.spatial.annotations import patch_spatial_annotations_from_table


def _make_table(instance_key="cell_id", ids=None):
    if ids is None:
        ids = [f"real_{index}" for index in range(5)]
    obs = pd.DataFrame(
        {
            instance_key: ids,
            "region": ["cell_circles"] * len(ids),
        },
        index=pd.Index([str(index) for index in range(len(ids))]),
    )
    adata = AnnData(
        X=np.ones((len(ids), 1), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(["gene"])),
    )
    adata.uns["spatialdata_attrs"] = {
        "region": "cell_circles",
        "region_key": "region",
        "instance_key": instance_key,
    }
    return adata


def _make_project(tmp_path, key_field="cell_id", key_values=None, extra_columns=None):
    if key_values is None:
        key_values = [f"{index}_0" for index in range(5)]
    data = {
        key_field: key_values,
        "value": list(range(len(key_values))),
    }
    if extra_columns:
        data.update(extra_columns)
    project_dir = tmp_path / "project"
    project = MDVProject(str(project_dir), delete_existing=True)
    project.add_datasource(
        "cells",
        pd.DataFrame(data),
        columns=[{"name": key_field, "field": key_field, "datatype": "unique"}],
    )
    return project_dir


def _annotation(ids=None):
    if ids is None:
        ids = ["real_0", "real_1", "real_2", "real_3", "real_4"]
    return pd.DataFrame(
        {
            "cell_id": ids,
            "celltype": ["A", "B", "C", "D", "E"][: len(ids)],
        }
    )


def test_patch_spatial_annotations_dry_run_does_not_write(tmp_path):
    project_dir = _make_project(tmp_path)
    adata = _make_table()

    report = patch_spatial_annotations_from_table(
        str(project_dir),
        _annotation(),
        adata,
        table_name="table",
    )

    project = MDVProject(str(project_dir))
    assert report.applied is False
    assert report.key_action == "repair"
    assert report.preservation_column == "cell_index"
    assert project.get_column("cells", "cell_id") == [f"{index}_0" for index in range(5)]
    fields = [column["field"] for column in project.get_datasource_metadata("cells")["columns"]]
    assert "cell_index" not in fields
    assert "celltype" not in fields


def test_patch_spatial_annotations_apply_preserves_old_key_and_adds_annotations(tmp_path):
    project_dir = _make_project(tmp_path)
    adata = _make_table()

    report = patch_spatial_annotations_from_table(
        str(project_dir),
        _annotation(ids=["real_0", "real_2", "real_4"]),
        adata,
        table_name="table",
        apply_changes=True,
    )

    project = MDVProject(str(project_dir))
    assert report.applied is True
    assert report.missing_annotation_count == 2
    assert report.backup_dir is not None
    assert os.path.exists(os.path.join(report.backup_dir, "datafile.h5"))
    assert os.path.exists(os.path.join(report.backup_dir, "datasources.json"))
    assert project.get_column("cells", "cell_id") == [f"real_{index}" for index in range(5)]
    assert project.get_column("cells", "cell_index") == [f"{index}_0" for index in range(5)]
    assert project.get_column("cells", "celltype") == ["A", "ND", "B", "ND", "C"]


def test_patch_spatial_annotations_is_idempotent(tmp_path):
    project_dir = _make_project(tmp_path)
    adata = _make_table()

    patch_spatial_annotations_from_table(
        str(project_dir),
        _annotation(),
        adata,
        table_name="table",
        apply_changes=True,
    )
    report = patch_spatial_annotations_from_table(
        str(project_dir),
        _annotation(),
        adata,
        table_name="table",
        apply_changes=True,
    )

    project = MDVProject(str(project_dir))
    fields = [column["field"] for column in project.get_datasource_metadata("cells")["columns"]]
    assert report.key_action == "unchanged"
    assert fields.count("cell_id") == 1
    assert fields.count("cell_index") == 1
    assert fields.count("celltype") == 1


def test_patch_spatial_annotations_rejects_duplicate_annotation_ids(tmp_path):
    project_dir = _make_project(tmp_path)
    adata = _make_table()
    annotations = pd.DataFrame(
        {
            "cell_id": ["real_0", "real_0"],
            "celltype": ["A", "B"],
        }
    )

    with pytest.raises(ValueError, match="duplicate"):
        patch_spatial_annotations_from_table(
            str(project_dir),
            annotations,
            adata,
            table_name="table",
        )


def test_patch_spatial_annotations_rejects_extra_annotation_ids(tmp_path):
    project_dir = _make_project(tmp_path)
    adata = _make_table()
    annotations = pd.DataFrame(
        {
            "cell_id": ["real_0", "not_in_project"],
            "celltype": ["A", "B"],
        }
    )

    with pytest.raises(ValueError, match="not present"):
        patch_spatial_annotations_from_table(
            str(project_dir),
            annotations,
            adata,
            table_name="table",
        )


def test_patch_spatial_annotations_rejects_preservation_conflict(tmp_path):
    project_dir = _make_project(
        tmp_path,
        extra_columns={"cell_index": [f"existing_{index}" for index in range(5)]},
    )
    adata = _make_table()

    with pytest.raises(ValueError, match="cell_index"):
        patch_spatial_annotations_from_table(
            str(project_dir),
            _annotation(),
            adata,
            table_name="table",
            apply_changes=True,
        )


def test_patch_spatial_annotations_supports_non_cell_id_instance_key(tmp_path):
    source_ids = [f"instance_{index}" for index in range(5)]
    project_dir = _make_project(
        tmp_path,
        key_field="instance_id",
        key_values=[f"old_{index}" for index in range(5)],
    )
    adata = _make_table(instance_key="instance_id", ids=source_ids)
    annotations = pd.DataFrame(
        {
            "instance_id": ["instance_0", "instance_1", "instance_2"],
            "celltype": ["A", "B", "C"],
        }
    )

    report = patch_spatial_annotations_from_table(
        str(project_dir),
        annotations,
        adata,
        table_name="table",
        apply_changes=True,
    )

    project = MDVProject(str(project_dir))
    assert report.instance_key == "instance_id"
    assert report.preservation_column == "instance_id_index"
    assert project.get_column("cells", "instance_id") == source_ids
    assert project.get_column("cells", "instance_id_index") == [f"old_{index}" for index in range(5)]
    assert project.get_column("cells", "celltype") == ["A", "B", "C", "ND", "ND"]


def test_patch_spatial_annotations_can_rename_incoming_annotation_key(tmp_path):
    project_dir = _make_project(tmp_path)
    adata = _make_table()
    annotations = pd.DataFrame(
        {
            "barcode": ["real_0", "real_1"],
            "celltype": ["A", "B"],
        }
    )

    patch_spatial_annotations_from_table(
        str(project_dir),
        annotations,
        adata,
        table_name="table",
        annotation_key="barcode",
        apply_changes=True,
    )

    project = MDVProject(str(project_dir))
    assert project.get_column("cells", "celltype") == ["A", "B", "ND", "ND", "ND"]
