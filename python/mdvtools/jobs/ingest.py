import h5py
import json


def ingest_column_output(project, params: dict, ws) -> dict | None:
    """Output shape column (ds, cols): write the worker's result column on to the datasouce
    set_column is idempotent - it deletes any existing column of the same name before adding, and updates-or-appends
    the metadata entry. So re-ingesting the same output_name REPLACES rather than duplicates (ADR0006).

    Returns:
        {"manifest": <worker manifest or None>, "outputs": [(datasource, column), ...]}.
    """
    out_name = params["output_name"]
    datasource = params["datasource"]
    with h5py.File(ws.output / "result.h5", "r") as f:
        ds = f[out_name]
        assert isinstance(ds, h5py.Dataset)
        values = [x.decode() for x in ds[:]]  # dataset comes back as bytes
    project.set_column(datasource, out_name, values)
    manifest = ws.output / "manifest.json"
    return {
        "manifest": json.loads(manifest.read_text()) if manifest.exists() else None,
        "outputs": [(datasource, out_name)],
    }
