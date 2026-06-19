from pathlib import Path
import h5py
import numpy as np

MARKER = "STATUS"  # worker writes "done/failed" - the primary completion signal


class Workspace:
    """Fixed per-job layout: input/ (tray), work/ (scratch), output/ (results).
    jobs_root is a sibling of any project export dir."""

    def __init__(self, jobs_root: Path, job_id: str):
        self.root = Path(jobs_root) / job_id
        self.input = self.root / "input"
        self.work = self.root / "work"
        self.output = self.root / "output"

        for d in (self.input, self.work, self.output):
            d.mkdir(parents=True, exist_ok=True)

    def read_marker(self) -> str | None:
        p = self.root / MARKER
        return p.read_text().strip() if p.exists() else None


# materialize keyed by INPUT shape, not tool
def materialize_columns_tray(project, spec, params: dict, ws: Workspace) -> None:
    """Input shape 'columns': every type == 'column' param -> a string dataset named by the param;
    scalar params (separator, output_name, ...) -> attrs. Shared by all column-input tools"""

    col_params = [p for p in spec.params if p.type == "column"]
    s = h5py.string_dtype()
    with project.lock("read"):
        cols = {
            p.name: project.get_column(params[p.options_from], params[p.name])
            for p in col_params
        }
    with h5py.File(ws.input / "tray.h5", "w") as f:
        for name, values in cols.items():
            f.create_dataset(name, data=np.array([str(x) for x in values], dtype=s))
        for p in spec.params:
            if p.type != "column" and params.get(p.name) is not None:
                f.attrs[p.name] = params[p.name]
