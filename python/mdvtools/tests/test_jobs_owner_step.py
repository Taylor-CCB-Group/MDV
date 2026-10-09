import json
import subprocess
import sys
from dataclasses import asdict

import h5py
import pandas as pd
from mdvtools.mdvproject import MDVProject
from mdvtools.jobs.jobstore import JobRecord


def _make_project(tmp_path):
    project = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    df = pd.DataFrame({"sample": ["s1", "s2", "s3"], "cluster": ["a", "b", "a"]})
    project.add_datasource("cells", df)
    return project


def _run_owner_step(step, project, ws_dir, rec):
    """Run an owner step the way the driver does (ADR-0014): its own process, record on stdin."""
    return subprocess.run(
        [sys.executable, "-m", "mdvtools.jobs.owner_step", step, project.dir, str(ws_dir)],
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
        assert [x.decode() for x in f["column_a"][:]] == ["s1", "s2", "s3"]
        assert [x.decode() for x in f["column_b"][:]] == ["a", "b", "a"]
        assert f.attrs["separator"] == "-"
