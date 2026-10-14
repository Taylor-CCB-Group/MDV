import json
import subprocess
import sys
from dataclasses import asdict

import h5py
import numpy as np
import pandas as pd
import scipy.sparse
from mdvtools.mdvproject import MDVProject
from mdvtools.jobs.jobstore import JobRecord


def _make_project(tmp_path):
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    df = pd.DataFrame({"sample": ["s1", "s2", "s3"], "cluster": ["a", "b", "a"]})
    project.add_datasource("cells", df)
    return project


def _run_owner_step(step, project, ws_dir, rec, rows_path=None):
    """Run an owner step the way the driver does (ADR-0014): its own process, record on stdin."""
    extra = [str(rows_path)] if rows_path is not None else []
    return subprocess.run(
        [sys.executable, "-m", "mdvtools.jobs.owner_step", step, project.dir, str(ws_dir), *extra],
        input=json.dumps(asdict(rec)),
        capture_output=True,
        text=True,
    )


def test_materialize_step_writes_the_tray_in_its_own_process(tmp_path):
    project = _make_project(tmp_path)
    rec = JobRecord(
        job_id="job1",
        tool_id="concat_columns",
        params={"datasource": "cells", "column_a": "sample", "column_b": "cluster",
                "separator": "-", "output_name": "out"},
    )
    ws_dir = tmp_path / "scratch" / rec.job_id

    result = _run_owner_step("materialize", project, ws_dir, rec)

    assert result.returncode == 0, result.stderr
    with h5py.File(ws_dir / "input" / "tray.h5", "r") as f:
        column_a, column_b = f["column_a"], f["column_b"]
        assert isinstance(column_a, h5py.Dataset) and isinstance(column_b, h5py.Dataset)
        assert [x.decode() for x in column_a[:]] == ["s1", "s2", "s3"]
        assert [x.decode() for x in column_b[:]] == ["a", "b", "a"]
        assert f.attrs["separator"] == "-"


def test_ingest_step_writes_columns_and_pointer_and_prints_provenance(tmp_path):
    project = _make_project(tmp_path)
    rec = JobRecord(
        job_id="job2",
        tool_id="concat_columns",
        params={"datasource": "cells", "column_a": "sample", "column_b": "cluster",
                "separator": "-", "output_name": "out"},
        created=100.0,
    )
    ws_dir = tmp_path / "scratch" / rec.job_id
    output = ws_dir / "output"
    output.mkdir(parents=True)
    with h5py.File(output / "result.h5", "w") as f:   # what the worker leaves behind
        f.create_dataset("out", data=["x", "y", "z"], dtype=h5py.string_dtype())
    (output / "manifest.json").write_text(json.dumps({"rows": 3}))

    result = _run_owner_step("ingest", project, ws_dir, rec)

    assert result.returncode == 0, result.stderr
    provenance = json.loads(result.stdout)
    assert provenance["job_id"] == "job2"
    assert provenance["submitted_at"] == 100.0
    assert provenance["output"] == {"rows": 3}
    assert project.get_column("cells", "out") == ["x", "y", "z"]
    assert project.get_column_metadata("cells", "out")["provenance"] == {
        "kind": "job",
        "job_id": "job2",
        "tool_id": "concat_columns",
        "content_hash": provenance["content_hash"],
    }


def test_materialize_step_cuts_the_rows_named_by_the_sidecar(tmp_path):
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    project.add_datasource("cells", pd.DataFrame({"cell_id": ["c0", "c1", "c2"]}))
    project.add_datasource("genes", pd.DataFrame({"name": ["g0", "g1"]}))
    project.add_rows_as_columns_link("cells", "genes", "name", "Gene Expr")
    X = scipy.sparse.csc_matrix(np.array([[1.0, 0.0], [0.0, 2.0], [3.0, 4.0]], dtype=np.float32))
    project.add_rows_as_columns_subgroup("cells", "genes", "gs", X, name="gene_scores")
    rows_path = tmp_path / "rows.npy"
    np.save(rows_path, np.array([1, 2], dtype=np.uint32))
    rec = JobRecord(
        job_id="job3",
        tool_id="umap",
        params={"datasource": "cells", "layer": "gs", "selection": "T cells", "output_name": "UMAP_T_cells"},
        selection_rows=2,
    )
    ws_dir = tmp_path / "scratch" / rec.job_id

    result = _run_owner_step("materialize", project, ws_dir, rec, rows_path)

    assert result.returncode == 0, result.stderr
    with h5py.File(ws_dir / "input" / "tray.h5", "r") as f:
        assert f.attrs["n_cells"] == 2


def test_ingest_step_spreads_a_selection_and_labels_the_pointer(tmp_path):
    project = _make_project(tmp_path)   # 3 cells
    rows_path = tmp_path / "rows.npy"
    np.save(rows_path, np.array([0, 2], dtype=np.uint32))
    rec = JobRecord(
        job_id="job4",
        tool_id="umap",
        params={"datasource": "cells", "layer": "gs", "selection": "T cells", "output_name": "UMAP_T_cells"},
        selection_rows=2,
    )
    ws_dir = tmp_path / "scratch" / rec.job_id
    output = ws_dir / "output"
    output.mkdir(parents=True)
    with h5py.File(output / "result.h5", "w") as f:
        f.create_dataset("UMAP_T_cells_1", data=np.array([0.5, 1.5]))

    result = _run_owner_step("ingest", project, ws_dir, rec, rows_path)

    assert result.returncode == 0, result.stderr
    values = np.asarray(project.get_column("cells", "UMAP_T_cells_1"), dtype=np.float64)
    assert values[0] == 0.5 and np.isnan(values[1]) and values[2] == 1.5
    pointer = project.get_column_metadata("cells", "UMAP_T_cells_1")["provenance"]
    assert pointer["selection"] == {"name": "T cells", "n_rows": 2}
