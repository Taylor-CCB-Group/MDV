from typing import cast

import h5py
import json
import numpy as np


def ingest_column_output(project, params: dict, ws, rows: np.ndarray | None = None) -> dict:
    """Output shape column (ds, cols): ingest every dataset the worker wrote to result.h5 as a column on the datasource.

    With `rows`, the sorted indices of a selection (ADR-0017), the worker's values fill those rows and
    every other row is NaN, MDV's missing value, so the column stays full length.

    Returns:
        {"manifest": <worker manifest or None>, "outputs": [(datasource, column), ...]}.
    """
    datasource = params["datasource"]
    size = int(project.get_datasource_metadata(datasource)["size"])
    outputs = []
    with h5py.File(ws.output / "result.h5", "r") as f:
        for col in sorted(f.keys()):
            ds = cast(h5py.Dataset, f[col])
            values = ds[:]
            if h5py.check_string_dtype(ds.dtype):
                values = [x.decode() for x in values] # string comes back as bytes -> text
            if rows is not None:
                full = np.full(size, np.nan)
                full[rows] = values
                values = full

            # numeric columns pass through as numpy array -> set_column infers double
            project.set_column(datasource, col, values)
            outputs.append((datasource, col))

    manifest = ws.output / "manifest.json"
    return {
        "manifest": json.loads(manifest.read_text()) if manifest.exists() else None,
        "outputs": outputs,
    }
