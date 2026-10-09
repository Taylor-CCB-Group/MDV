# using the courier model described in ADR-0004
from pathlib import Path
import h5py
import numpy as np
import json
import traceback


def run(workspace: str) -> None:
    ws = Path(workspace)
    try:
        with h5py.File(ws / "input" / "tray.h5", "r") as f:
            col_a = f["column_a"]
            col_b = f["column_b"]
            assert isinstance(col_a, h5py.Dataset)
            assert isinstance(col_b, h5py.Dataset)
            a = [x.decode() for x in col_a[:]]  # dataset string comes back as bytes
            b = [x.decode() for x in col_b[:]]
            sep = f.attrs["separator"]
            out_name = f.attrs["output_name"]
        joined = [f"{x}{sep}{y}" for x, y in zip(a, b)]

        with h5py.File(ws / "output" / "result.h5", "w") as f:
            f.create_dataset(out_name, data=np.array(joined, dtype=h5py.string_dtype()))
            f.attrs["output_name"] = out_name
        (ws / "output" / "manifest.json").write_text(
            json.dumps({"rows": len(joined)})
        )  # basic provenance
        (ws / "STATUS").write_text("done")
    except Exception as e:
        (ws / "STATUS").write_text("failed")
        (ws / "output" / "error.txt").write_text(traceback.format_exc())
        raise
