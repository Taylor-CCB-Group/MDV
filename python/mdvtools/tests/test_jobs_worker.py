import subprocess
import sys

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
        ds = f["joined"]
        assert isinstance(ds, h5py.Dataset)
        values = [x.decode() for x in ds[:]]
    assert values == ["a-1", "b-2", "c-3"]
    assert (ws.output / "manifest.json").exists()


def test_worker_writes_failed_marker_when_tray_missing(tmp_path):
    ws = Workspace(tmp_path, "job2")  # no tray.h5 built
    with pytest.raises(Exception):
        run(str(ws.root))
    assert ws.read_marker() == "failed"
    assert (ws.output / "error.txt").exists()


def test_worker_import_does_not_load_mdvproject():
    # A compute node runs `python -m mdvtools.jobs.run_worker <module:fn> <ws>`, which imports
    # the `mdvtools` package first. A fresh interpreter shows whether that drags in MDVProject
    # and its stack (pandas, scanpy) when the worker itself only needs h5py and numpy.
    code = (
        "import sys\n"
        "import mdvtools.jobs.run_worker\n"
        "import mdvtools.jobs.workers.concat_worker\n"
        "print(sorted(m for m in ('mdvtools.mdvproject', 'pandas') if m in sys.modules))\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    ).stdout
    assert out.strip() == "[]"


def test_top_level_exports_still_resolve():
    import mdvtools
    from mdvtools.mdvproject import MDVProject

    assert mdvtools.MDVProject is MDVProject
    assert callable(mdvtools.setup_ligand_network)
