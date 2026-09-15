"""Demo / inspection for the jobs framework — pick a job (concat_columns or umap) and watch it run
end-to-end: durable record first, tray staged, worker subprocess, ingest, provenance pointer."""

import json
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse

from mdvtools.mdvproject import MDVProject
from mdvtools.jobs.manager import JobManager


def rule(title):
    print("\n" + "=" * 72 + f"\n{title}\n" + "=" * 72)


def tree(root):
    if not root.exists():
        print(f"  {root}  ->  (gone)")
        return
    for p in sorted(root.rglob("*")):
        depth = len(p.relative_to(root).parts) - 1
        print(f"  {'    ' * depth}{p.name}{'/' if p.is_dir() else ''}")


def dump_h5(path):
    if not path.exists():
        print(f"    (no {path.name})")
        return
    with h5py.File(path, "r") as f:
        for k in f.keys():
            vals = [x.decode() if isinstance(x, bytes) else x for x in f[k][:]]
            if len(vals) > 12:                       # keep the CSC triplet readable
                vals = vals[:6] + [f"... ({len(vals)} total)"]
            print(f"    dataset {k!r}: {vals}")
        for k, v in f.attrs.items():
            print(f"    attr    {k!r}: {v!r}")


def show_json(obj):
    print("    " + json.dumps(obj, indent=2).replace("\n", "\n    "))


def ask(label, default, cast):
    raw = input(f"    {label} [{default}]: ").strip()
    return cast(raw) if raw else default


def build_concat(project):
    """Columns-input job: join two text columns on 'cells'."""
    project.add_datasource(
        "cells", pd.DataFrame({"sample": ["s1", "s2", "s3"], "cluster": ["a", "b", "a"]})
    )
    return "concat_columns", {
        "datasource": "cells",
        "column_a": "sample",
        "column_b": "cluster",
        "separator": "_",
        "output_name": "sample_cluster",
    }


def build_umap(project):
    """Matrix-input job: a stored sparse gs expression matrix; prompt the scanpy knobs."""
    n_cells, n_genes = 60, 10
    project.add_datasource("cells", pd.DataFrame({"cell_id": [f"c{i}" for i in range(n_cells)]}))
    project.add_datasource("genes", pd.DataFrame({"name": [f"g{j}" for j in range(n_genes)]}))
    project.add_rows_as_columns_link("cells", "genes", "name", "Gene Expr")
    rng = np.random.default_rng(0)
    dense = rng.random((n_cells, n_genes)).astype(np.float32)
    dense[dense < 0.6] = 0.0                                    # genuinely sparse
    project.add_rows_as_columns_subgroup(
        "cells", "genes", "gs", scipy.sparse.csc_matrix(dense), name="gene_scores"
    )
    print("\n  UMAP parameters (Enter for the scanpy default):")
    return "umap", {
        "datasource": "cells",
        "layer": "gs",
        "output_name": "UMAP",
        "n_neighbors": ask("n_neighbors", 15, int),
        "min_dist": ask("min_dist", 0.5, float),
        "n_components": ask("n_components", 2, int),
    }


def main():
    base = Path(tempfile.mkdtemp(prefix="mdv_jobs_demo_"))
    project = MDVProject(str(base / "proj"), delete_existing=True)
    scratch_root = base / "scratch"      # stands in for HPC/cluster scratch — OUTSIDE the project

    print("Which job?  [1] concat_columns   [2] umap")
    tool_id, params = build_umap(project) if input("> ").strip() == "2" else build_concat(project)
    datasource = params["datasource"]

    mgr = JobManager(project, workspace_root=scratch_root, max_concurrent=2)

    rule("0. Where everything lives  (two SEPARATE roots, linked only by job_id)")
    print(f"  project       : {project.dir}")
    print(f"  job records   : {mgr.records_root}/records   (durable, INSIDE the project)")
    print(f"  job workspaces: {mgr.workspace_root}   (ephemeral scratch, OUTSIDE the project)")

    rule("1. submit()  ->  durable record written FIRST (write-ahead), then staged + launched")
    before = {c["field"] for c in project.get_datasource_metadata(datasource)["columns"]}
    job_id = mgr.submit(tool_id, params)
    ws = mgr.workspace_root / job_id
    rec_file = mgr.records_root / "records" / f"{job_id}.json"
    print(f"  tool  : {tool_id}")
    print(f"  params: {params}")
    print(f"  job_id: {job_id}")
    print("\n  records/<job_id>.json  (status RUNNING, provenance still null):")
    show_json(json.loads(rec_file.read_text()))
    print("\n  input/tray.h5  — the owner staged this (decoded at the boundary; worker imports no MDV):")
    dump_h5(ws / "input" / "tray.h5")

    rule("2. wait for the worker subprocess  (owner has NOT touched the project yet)")
    marker = ws / "STATUS"
    deadline = time.time() + 240         # umap runs real scanpy in a subprocess (import/JIT on first run)
    while time.time() < deadline and not marker.exists():
        time.sleep(0.1)
    status = marker.read_text().strip() if marker.exists() else "(timed out)"
    print(f"  STATUS marker (written LAST): {status!r}")
    if status != "done":
        err = ws / "output" / "error.txt"
        print("\n  worker did not succeed — error.txt:")
        print("    " + (err.read_text() if err.exists() else "(none)").replace("\n", "\n    "))
        print(f"\n  Left for inspection: {base}")
        return
    print("\n  output/result.h5:")
    dump_h5(ws / "output" / "result.h5")
    print("\n  output/manifest.json:")
    print("    " + (ws / "output" / "manifest.json").read_text())

    rule("3. tick()  ->  ingest, promote provenance into the record, clean the workspace")
    mgr.tick()
    print(f"  workspace after success: {'GONE — cleaned (ADR-0007)' if not ws.exists() else 'STILL THERE'}")
    rec = json.loads(rec_file.read_text())
    print(f"  record status: {rec['status']!r};  provenance now lives in the record:")
    show_json(rec["provenance"])

    rule("4. the output column(s), each carrying a provenance POINTER (ADR-0009)")
    after = {c["field"] for c in project.get_datasource_metadata(datasource)["columns"]}
    for col in sorted(after - before):
        meta = next(
            c for c in project.get_datasource_metadata(datasource)["columns"] if c["field"] == col
        )
        pointer = meta.get("provenance")
        resolved = project.get_column_provenance(datasource, col)
        values = list(project.get_column(datasource, col))
        preview = values[:6] + ["..."] if len(values) > 6 else values
        print(f"\n  {col} ({meta['datatype']}): {preview}")
        print(f"    pointer on column : {pointer}")
        print(f"    resolves to record: job_id={resolved['job_id']}  _resolved={resolved['_resolved']}")
        print(
            f"    one hash, not two : column={pointer['content_hash']!r}  "
            f"record={resolved['provenance']['content_hash']!r}"
        )

    print(f"\n  Left for inspection: {base}")


if __name__ == "__main__":
    main()
