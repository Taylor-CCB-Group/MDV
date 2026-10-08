# Project JSON scripts

Scan column names across MDV projects, and snapshot the JSON files those scans read. Both tools walk directories that contain `datasources.json`. They do not open HDF5 or Zarr, do not touch the database, and do not import `mdvtools`. A directory with `datasources.json` is a project; the walk stops there, so `history/` and other nested data are not treated as further projects.

Run them from the repo root:

```bash
python python/mdvtools/scripts/column_name_report.py --help
python python/mdvtools/scripts/project_json_history.py --help
```

## Column name report

`column_name_report.py` groups each datasource's `field` values by spelling. Trimming, lowercasing, and dropping every character outside `[a-z0-9]` makes `cell_type`, `Cell Type`, and `cell-type` one cluster (`celltype`). Every cluster is reported, including a field that uses one spelling everywhere. A column is marked used when a chart in that project's `views.json` references its exact `field` on the same datasource. The report records the view name and chart type.

Markdown goes to stdout. `--csv` and `--json` also write those formats. CSV columns are `datasource`, `normalized`, `field`, `display_name`, `project`, `used_in_views`, `view_name`, and `chart_type`.

```bash
python python/mdvtools/scripts/column_name_report.py /path/to/projects
python python/mdvtools/scripts/column_name_report.py /path/to/projects --csv report.csv --json report.json
```

`--fuzzy RATIO` (0–1) also merges normalized keys in the same datasource whose `difflib` ratio is at least `RATIO`. It is off by default. Those matches are suggestions and can be wrong. Keys shorter than `--min-fuzzy-length` (default 6) stay separate, so `pc1` and `pc2` do not merge. Field names that start with `__` are omitted unless you pass `--include-internal`.

## JSON history

`project_json_history.py` copies `datasources.json` and, when present, `views.json` and `state.json` into `history/<revision>/` inside each project. Take a snapshot before editing those files, then restore a revision if the edit should be undone. Data files are left as they are, so a restored JSON file can refer to columns or views the data no longer matches.

One `snapshot` invocation stamps the same run id on every project. Revision ids differ per project. `log` lists revisions newest first, with the parent revision, run id, time, and message. Restore one project by revision id, or every project from one snapshot with `--run`.

```bash
python python/mdvtools/scripts/project_json_history.py snapshot /path/to/projects --message "before rename"
python python/mdvtools/scripts/project_json_history.py log /path/to/projects
python python/mdvtools/scripts/project_json_history.py restore /path/to/one/project <revision-id>
python python/mdvtools/scripts/project_json_history.py restore /path/to/projects --run <run-id>
```

`snapshot` prints the run id, then each project path and its new revision id.

## Tests

```bash
./venv/bin/python -m pytest python/mdvtools/tests/test_column_name_report.py python/mdvtools/tests/test_project_json_history.py
```

## Other scripts

This directory also contains `migrate_projects_autoincrement.py` (rebuild a SQLite database so project ids are not reused) and `manage_project_permissions.py` (grant project access). Each script's `--help` and module docstring describe that command.
