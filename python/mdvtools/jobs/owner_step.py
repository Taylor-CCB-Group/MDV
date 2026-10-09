"""Owner steps (ADR-0014): the owner-side work that reads or writes the project store for one job.

The driver runs each step in its own short-lived process and waits on it, so the HDF5 work never
runs on the server's gevent hub or beside its request handlers in the same process.

Usage: python -m mdvtools.jobs.owner_step <materialize|ingest> <project_dir> <workspace_dir>
The job record arrives as JSON on stdin.
"""
import json
import sys
from pathlib import Path

from mdvtools.mdvproject import MDVProject
from .jobstore import JobRecord
from .registry import get_tool
from .workspace import Workspace, materialize_columns_tray, materialize_matrix_tray

# materializers keyed by INPUT shape
MATERIALIZERS = {"columns": materialize_columns_tray, "matrix": materialize_matrix_tray}


def materialize(project, rec: JobRecord, ws: Workspace) -> None:
    spec = get_tool(rec.tool_id)
    MATERIALIZERS[spec.input_shape](project, spec, rec.params, ws)


def main(argv: list[str]) -> None:
    step, project_dir, workspace_dir = argv
    rec = JobRecord(**json.loads(sys.stdin.read()))
    project = MDVProject(project_dir)
    ws_path = Path(workspace_dir)
    ws = Workspace(ws_path.parent, ws_path.name)
    if step == "materialize":
        materialize(project, rec, ws)
    else:
        raise ValueError(f"unknown owner step {step!r}")


if __name__ == "__main__":
    main(sys.argv[1:])
