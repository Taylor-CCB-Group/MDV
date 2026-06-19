import pytest
import h5py
import numpy as np

from mdvtools.jobs.workspace import Workspace
from mdvtools.jobs.workers.concat_worker import run


def _build_tray(ws, col_a, col_b, separator, output_name):
    s = h5py.string_dtype()
    with h5py.File(ws.input / "tray.h5", "w") as f:
        f.create_dataset("column_a", data=np.array(col_a, dtype=s))
        f.create_dataset("column_b", data=np.array(col_b, dtype=s))
        f.attrs["separator"] = separator
        f.attrs["output_name"] = output_name


def test_worker_joins_column(tmp_path):
    # No MDVProject anywhere - proves the worker is isolated from mdvtools (courier model)
    ws = Workspace(tmp_path, "job1")
    _build_tray(ws, ["a", "b", "c"], ["1", "2", "3"], "-", "joined")

    run(str(ws.root))

    assert ws.read_marker() == "done"
    with h5py.File(ws.output / "result.h5", "r") as f:
        assert f.attrs["output_name"] == "joined"
        values = [x.decode() for x in f["joined"][:]]
    assert values == ["a-1", "b-2", "c-3"]
    assert (ws.output / "manifest.json").exists()


def test_worker_writes_failed_marker_when_tray_missing(tmp_path):
    ws = Workspace(tmp_path, "job2")  # no tray.h5 built
    with pytest.raises(Exception):
        run(str(ws.root))
    assert ws.read_marker() == "failed"
    assert (ws.output / "error.txt").exists()
