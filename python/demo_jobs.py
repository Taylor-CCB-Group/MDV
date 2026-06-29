"""Demo / inspection for the jobs framework (concat_columns)."""
import json
import tempfile
import time
from pathlib import Path

import h5py
import pandas as pd

from mdvtools.mdvproject import MDVProject
from mdvtools.jobs.manager import JobManager


def rule(title: str) -> None:
    print("\n" + "=" * 72 + f"\n{title}\n" + "=" * 72)


def tree(root: Path) -> None:
    if not root.exists():
        print(f"  {root}  ->  (gone)")
        return
    for p in sorted(root.rglob("*")):
        depth = len(p.relative_to(root).parts) - 1
        print(f"  {'    ' * depth}{p.name}{'/' if p.is_dir() else ''}")


def dump_h5(path: Path) -> None:
    if not path.exists():
        print(f"    (no {path.name})")
        return
    with h5py.File(path, "r") as f:
        for k in f.keys():
            vals = [x.decode() if isinstance(x, bytes) else x for x in f[k][:]]
            print(f"    dataset {k!r}: {vals}")
        for k, v in f.attrs.items():
            print(f"    attr    {k!r}: {v!r}")


def show_json(obj) -> None:
    print("    " + json.dumps(obj, indent=2).replace("\n", "\n    "))


def main() -> None:
    base = Path(tempfile.mkdtemp(prefix="mdv_jobs_demo_"))
    proj_dir = base / "proj"
    scratch_root = base / "scratch"  # stands in for HPC/cluster scratch — OUTSIDE the project

    project = MDVProject(str(proj_dir), delete_existing=True)
    project.add_datasource("cells", pd.DataFrame({
        "sample": ["s1", "s2", "s3"],
        "cluster": ["a", "b", "a"],
    }))
    # records_root defaults to <project>/jobs (inside the project); workspace_root is the
    # ephemeral scratch we point outside the project, as HPC would.
    mgr = JobManager(project, workspace_root=scratch_root, max_concurrent=2)

    rule("0. Where everything lives  (two SEPARATE roots, linked only by job_id)")
    print(f"  demo base     : {base}")
    print(f"  project       : {proj_dir}   (datafile.h5 + datasources.json)")
    print(f"  job records   : {mgr.records_root}/records   (durable, INSIDE the project — travels with it,")
    print( "                  the catalog never scans a project's subdir, exporter skips it)")
    print(f"  job workspaces: {mgr.workspace_root}   (ephemeral scratch, OUTSIDE the project — on HPC: $SCRATCH)")
    print( "  link          : records/<job_id>.json  <->  <scratch>/<job_id>/   (same job_id, two locations)")

    rule("1. submit()  ->  durable record written FIRST (write-ahead), then staged + launched")
    job_id = mgr.submit("concat_columns", {
        "datasource": "cells",
        "column_a": "sample",
        "column_b": "cluster",
        "separator": "_",
        "output_name": "sample_cluster",
    })
    ws = mgr.workspace_root / job_id
    rec_file = mgr.records_root / "records" / f"{job_id}.json"
    print(f"  job_id: {job_id}")
    print("\n  records/<job_id>.json  (status RUNNING, handle set, provenance still null):")
    show_json(json.loads(rec_file.read_text()))
    print(f"\n  workspace {ws}:")
    tree(ws)
    print("\n  input/tray.h5  — the owner staged this (columns already decoded to real values):")
    dump_h5(ws / "input" / "tray.h5")

    rule("2. wait for the worker subprocess  (owner has NOT touched the project yet)")
    marker = ws / "STATUS"
    deadline = time.time() + 30
    while time.time() < deadline and not marker.exists():
        time.sleep(0.05)
    print(f"  STATUS marker (written LAST): {marker.read_text().strip()!r}")
    print("\n  output/ now holds the worker's answer:")
    tree(ws / "output")
    print("\n  output/result.h5:")
    dump_h5(ws / "output" / "result.h5")
    print("\n  output/manifest.json  (the provenance seed the worker emits):")
    print("    " + (ws / "output" / "manifest.json").read_text())

    mgr.tick()
    print(f"  workspace after success: {'GONE — cleaned (ADR-0007)' if not ws.exists() else 'STILL THERE'}")
    rec = json.loads(rec_file.read_text())
    print(f"\n  records/<job_id>.json  (status {rec['status']!r}).")
    print("  >>> PROVENANCE LIVES HERE — in the JSON record, not a .txt file:")
    show_json(rec["provenance"])

    rule("3. tick()  ->  ingest, promote provenance into the record, clean the workspace")
    rule("4. the output column, now part of the project")
    print(f"  get_column('cells', 'sample_cluster') -> {project.get_column('cells', 'sample_cluster')}")
    col = next(c for c in project.get_datasource_metadata("cells")["columns"]
               if c["field"] == "sample_cluster")
    print("\n  its entry in datasources.json (note: NO provenance on the column yet — that's the")
    print("  open discussion item; provenance currently lives ONLY in the job record above):")
    show_json(col)

    print(f"\n  Left for inspection: {base}")
    print(f"   - {mgr.records_root}/records/{job_id}.json   (record + provenance, durable, INSIDE the project)")
    print(f"   - {proj_dir}/datasources.json          (the new column)")
    print(f"   - {mgr.workspace_root}/{job_id}   is GONE (scratch cleaned on success); its contents shown above.")

if __name__ == "__main__":
    main()
