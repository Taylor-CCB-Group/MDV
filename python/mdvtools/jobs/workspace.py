from pathlib import Path
import tempfile
import h5py
import numpy as np
import json

MARKER = "STATUS"  # worker writes "done/failed" - the primary completion signal


def default_workspace_root(project) -> Path:
    """Default ephemeral-scratch root for per-job workspaces (ADR-0007).

    Lives OUTSIDE the project — never under the catalog's projects base_dir — so neither
    the catalog scan nor the project exporter ever sees it. Namespaced by project so two
    projects can't collide. Override on HPC to point at `$SCRATCH` / a shared mount the
    compute node can reach. The per-job dir under it is keyed by job_id, the only link
    back to the durable owner-side record (which lives inside the project)."""
    return Path(tempfile.gettempdir()) / "mdv_job_workspaces" / Path(project.dir).name


class Workspace:
    """Fixed per-job layout under the ephemeral scratch root (ADR-0007): input/ (tray),
    work/ (worker scratch), output/ (results). The scratch root lives OUTSIDE the project
    — on HPC it is cluster scratch the worker reaches. The per-job dir is keyed by job_id,
    the sole link to the durable owner-side record."""

    def __init__(self, workspace_root: Path, job_id: str):
        self.root = Path(workspace_root) / job_id
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

def materialize_matrix_tray(project, spec, params: dict, ws: Workspace) -> None:
    """Input shape 'matrix': copy the expression matrix (a rows_as_column subgroup) into its tray
    as its CSC triplet - x/i/p unchanged, plus explicit shape, stays sparse. The worker rebuilds csc_matrix
    from the triplet and imports no MDV.
    """
    datasource = params["datasource"]
    subgroup_key = params.get("layer", "gs")
    name, is_sparse = project._resolve_rows_as_columns_subgroup(datasource, subgroup_key)
    if not is_sparse:
        raise NotImplementedError("matrix tray materialization only supports sparse matrices")

    n_cells = int(project.get_datasource_metadata(datasource)["size"])
    with project.lock("read"):
        h5 = project._get_h5_handle(read_only = True)
        try:
            grp = h5[f"{datasource}/{name}"]
            assert isinstance(grp, h5py.Group)
            x = np.asarray(grp["x"], dtype=np.float32)
            i = np.asarray(grp["i"], dtype=np.uint32)
            p = np.asarray(grp["p"])
        finally:
            h5.close()
    n_genes = len(p) - 1

    with h5py.File(ws.input / "tray.h5", "w") as f:
        f.create_dataset("x", data=x)
        f.create_dataset("i", data=i)
        f.create_dataset("p", data=p)
        f.attrs["n_cells"] = n_cells
        f.attrs["n_genes"] = n_genes
        f.attrs["sparse"] = True
        # route params by their declared target call: control/output params (applies_to = None) stay
        # flat attrs; compute params bundle into kwargs.<call>, so the worker splats them into the call
        kwargs_bundles: dict[str, dict] = {}
        for pspec in spec.params:
            if pspec.type == "column" or params.get(pspec.name) is None:
                continue
            if pspec.applies_to is None:
                f.attrs[pspec.name] = params[pspec.name]
            else:
                kwargs_bundles.setdefault(pspec.applies_to, {})[pspec.name] = params[pspec.name]
        for group, kw in kwargs_bundles.items():
            f.attrs[f"kwargs.{group}"] = json.dumps(kw)
