# Curate an MDV project from the Data Cleaning sheet

`python/mdvtools/examples/curate_mdv_project.py` renames column display names, soft-deletes columns the sheet does not keep, and adds question views. Chart `param` values stay as field ids. `datafile.h5` is never modified.

Run it from the repository root with the project virtualenv. The workbook must contain a sheet named `Data Cleaning`. `--project-id` is the value in the sheet's `Project id` column. It can differ from the folder name.

The default is a dry run. It writes one JSON report beside the intended output directory and does not create that directory. `--write` copies the input project to `--output`, backs up `views.json`, and applies the plan there.

## Arguments

- `--project` — input MDV project directory (`datasources.json`, `views.json`, `state.json`, `datafile.h5`)
- `--project-id` — Excel `Project id`
- `--xlsx` — workbook path
- `--output` — directory for the curated project; must be a different directory from the input
- `--write` — apply the plan; omit it for a dry run

## What the dry run prints

- Planned renames, from `Final Column Name`. If Final is blank, `new_display_name` is used.
- If either of those names is repeated for the same project and datasource, that row uses `displayed name` instead. A blank or colliding displayed name keeps the column.
- Soft-deletes for columns the sheet does not keep. A column stays when a saved view or a new question view references it, the target is ambiguous, or the column is spatial. Soft-delete sets `deleted: true` and leaves the h5 data in place.
- Question views that would be added, and questions that are skipped. Protein concordance is skipped when there is no protein table. Spatial charts, collection-score columns, and DEG columns are always skipped. Collection charts plot member genes that already exist.
- The paths a later `--write` would use for `views.json.bak.<timestamp>` inside the output project and for the JSON backup beside it.

The report file is `<output directory name>_dry_run_report.json` in the parent of `--output`. After `--write` the file is named `<output directory name>_report.json`.

## Local project

Dry run for the original synovial project. The sheet id is `mdv_synovial`; the folder name is different.

```bash
./venv/bin/python python/mdvtools/examples/curate_mdv_project.py \
  --project ~/mdv/mdv_synovial_orig \
  --project-id mdv_synovial \
  --xlsx ~/mdv/MIM_Nov25.xlsx \
  --output /tmp/mdv_synovial_orig_curated
```

The report is `/tmp/mdv_synovial_orig_curated_dry_run_report.json`. `/tmp/mdv_synovial_orig_curated` is not created.

Apply the same plan into a new local folder:

```bash
./venv/bin/python python/mdvtools/examples/curate_mdv_project.py \
  --project ~/mdv/mdv_synovial_orig \
  --project-id mdv_synovial \
  --xlsx ~/mdv/MIM_Nov25.xlsx \
  --output ~/mdv/mdv_synovial_orig_curated \
  --write
```

`--write` copies `datafile.h5` into the output project and then checks that its size and modification time are unchanged. The input project's h5 is only opened for reading.

## Project on bmrc2

`bmrc2` in `~/.ssh/config` is `ucj494@cluster2.bmrc.ox.ac.uk`. Log in from a terminal that the cluster accepts, with this repository available on bmrc2 and a Python environment that can import `mdvtools`, pandas, and h5py. Reading the workbook also needs openpyxl.

Dry run for StromalAtlas. The output folder is not created.

```bash
ssh bmrc2 'python python/mdvtools/examples/curate_mdv_project.py \
  --project /gpfs3/well/cartography/shared/MDV/ondemand/projects/StromalAtlas_H52MDV \
  --project-id StromalAtlas_H52MDV \
  --xlsx /well/cartography/shared/MDV/MIM_Nov25.xlsx \
  --output /well/cartography/shared/MDV/StromalAtlas_H52MDV_curated'
```

Run that command from the repository root on the cluster, or pass the full path to `curate_mdv_project.py`. The report is written next to the output path:

`/well/cartography/shared/MDV/StromalAtlas_H52MDV_curated_dry_run_report.json`

Apply it after the report looks right:

```bash
ssh bmrc2 'python python/mdvtools/examples/curate_mdv_project.py \
  --project /gpfs3/well/cartography/shared/MDV/ondemand/projects/StromalAtlas_H52MDV \
  --project-id StromalAtlas_H52MDV \
  --xlsx /well/cartography/shared/MDV/MIM_Nov25.xlsx \
  --output /well/cartography/shared/MDV/StromalAtlas_H52MDV_curated \
  --write'
```

If `--output` already exists, `--write` stops and leaves it alone.
