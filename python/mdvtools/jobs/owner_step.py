"""Owner steps (ADR-0014): the owner-side work that reads or writes the project store for one job.

The driver runs each step in its own short-lived process and waits on it, so the HDF5 work never
runs on the server's gevent hub or beside its request handlers in the same process.

Usage: python -m mdvtools.jobs.owner_step <materialize|ingest> <project_dir> <workspace_dir> [rows]
The job record arrives as JSON on stdin. ingest prints the job's provenance as JSON on stdout.
`rows` is the path of the job's selection sidecar (ADR-0017), given only when it runs on a selection.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

from mdvtools.mdvproject import MDVProject
from .jobstore import JobRecord
from .registry import get_tool
from .workspace import Workspace, materialize_columns_tray, materialize_matrix_tray
from .ingest import ingest_column_output
from .provenance import build_provenance

# materializers keyed by INPUT shape; ingest handlers keyed by OUTPUT shape
MATERIALIZERS = {"columns": materialize_columns_tray, "matrix": materialize_matrix_tray}
INGESTERS = {"column": ingest_column_output}


def materialize(project, rec: JobRecord, ws: Workspace, rows: np.ndarray | None = None) -> None:
    spec = get_tool(rec.tool_id)
    MATERIALIZERS[spec.input_shape](project, spec, rec.params, ws, rows)


def ingest(project, rec: JobRecord, ws: Workspace) -> dict:
    """Write the worker's outputs and their provenance pointers into the project; return the provenance."""
    spec = get_tool(rec.tool_id)
    with project.lock("write"):
        result = INGESTERS[spec.output.shape](project, rec.params, ws)
        provenance = build_provenance(rec, result["manifest"])
        pointer = {
            "kind": "job",
            "job_id": rec.job_id,
            "tool_id": rec.tool_id,
            "content_hash": provenance["content_hash"],
        }
        for ds, col in result["outputs"]:
            project.set_column_metadata(ds, col, "provenance", pointer)
    return provenance


def main(argv: list[str]) -> None:
    step, project_dir, workspace_dir, *rows_path = argv
    # MDV's log handler writes to stdout, so stdout would mix log lines with the result. Point
    # fd 1 at stderr for the step and write the result to a copy of the original fd 1.
    sys.stdout.flush()
    result_out = os.fdopen(os.dup(1), "w")
    os.dup2(2, 1)

    rec = JobRecord(**json.loads(sys.stdin.read()))
    project = MDVProject(project_dir)
    ws_path = Path(workspace_dir)
    ws = Workspace(ws_path.parent, ws_path.name)
    rows = np.load(rows_path[0]) if rows_path else None
    if step == "materialize":
        materialize(project, rec, ws, rows)
    elif step == "ingest":
        result_out.write(json.dumps(ingest(project, rec, ws)))
    else:
        raise ValueError(f"unknown owner step {step!r}")
    result_out.close()


if __name__ == "__main__":
    main(sys.argv[1:])
