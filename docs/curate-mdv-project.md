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

## Project on a cluster

Log in to the cluster, with this repository available there and a Python environment that can import `mdvtools`, pandas, and h5py. Reading the workbook also needs openpyxl. Keep real hosts and paths out of this repository; the commands below use placeholders.

Dry run. The output folder is not created.

```bash
python python/mdvtools/examples/curate_mdv_project.py \
  --project <projects_dir>/<project> \
  --project-id <sheet_project_id> \
  --xlsx <workbook.xlsx> \
  --output <work_dir>/<project>_curated
```

Run that command from the repository root on the cluster, or pass the full path to `curate_mdv_project.py`. The report is written next to the output path:

`<work_dir>/<project>_curated_dry_run_report.json`

Apply it after the report looks right:

```bash
python python/mdvtools/examples/curate_mdv_project.py \
  --project <projects_dir>/<project> \
  --project-id <sheet_project_id> \
  --xlsx <workbook.xlsx> \
  --output <work_dir>/<project>_curated \
  --write
```

If `--output` already exists, `--write` stops and leaves it alone.

## Add the overview view

`python/mdvtools/template_views/create_overview_view.py` adds a `Cell atlas` view with the same charts and chart sizes as the one in `myeloid_v2_testing_latest`: the project title and abstract across the top, UMAPs coloured by cell type, cell group, tissue and disease, and composition charts below. It then orders the views: `Cell atlas` first, the question views in question order, then every other listed view in its current order. `initial_view` is set to `Cell atlas`. Only `views.json` and `state.json` change. `datafile.h5` is not opened.

Columns are chosen from the cells datasource by display name, after curation:

- cell type: `Cell State`, then `Cell Type`
- cell group: `Cell Subgroup`, then `Cell Group`
- disease: `Diagnosis`, then `Disease`
- tissue: `Tissue Simple`, then `Tissue`, then `Biopsy Site`
- UMAP: the field ending in `X_umap_mindist_0.25_1` / `_2`, then display names `UMAP 1` / `UMAP 2`

Deleted columns are skipped. If a role is not found, the script stops and names the flag to pass (`--cell-type`, `--cell-group`, `--disease`, `--tissue`, `--umap-x`, `--umap-y`). Each flag takes a field id or a display name.

Take a snapshot first, then dry run:

```bash
./venv/bin/python python/mdvtools/scripts/project_json_history.py snapshot <projects_dir>/mdv_synovial --message "before overview view"
./venv/bin/python python/mdvtools/template_views/create_overview_view.py \
  --project <projects_dir>/mdv_synovial
```

Put the abstract in a text file (Markdown is allowed), then write:

```bash
./venv/bin/python python/mdvtools/template_views/create_overview_view.py \
  --project <projects_dir>/mdv_synovial \
  --title "Synovial tissue atlas" \
  --abstract-file ~/synovial_abstract.md \
  --write
```

`--write` copies `views.json` and `state.json` to `*.bak.<timestamp>` beside them first. If `Cell atlas` already exists the script stops; pass `--replace` to rebuild it. Other keys in `state.json` are left as they are, including `show_gallery_on_open`, which opens the project on the view gallery instead of `initial_view` when true.

## Compound dot plots in question views

Where a question view compares groups (Q6 and Q7 by diagnosis, Q8 and Q9 by response or treatment stage), the scripts now add one `compound_dot_plot` split by that column, with a panel per group on one colour and size scale. They used to add a separate `dot_plot` per group, each restricted with a background filter. Single filtered plots, such as the Q2 Patient plot and the Q3 RNA and protein plots, stay `dot_plot`. Views with compound dot plots need an MDV build that includes the compound dot plot (main from 28 October 2026 onwards). Older builds can't open them.

For a project that already has the separate plots, `python/mdvtools/scripts/convert_to_compound_dot_plots.py` merges them in `views.json` and updates the "About this page" text to match. It changes nothing else. Take a snapshot, dry run, then write:

```bash
./venv/bin/python python/mdvtools/scripts/project_json_history.py snapshot <projects_dir>/mdv_synovial --message "before compound dot plots"
./venv/bin/python python/mdvtools/scripts/convert_to_compound_dot_plots.py \
  --project <projects_dir>/mdv_synovial
./venv/bin/python python/mdvtools/scripts/convert_to_compound_dot_plots.py \
  --project <projects_dir>/mdv_synovial --write
```

`--write` copies `views.json` to `views.json.bak.<timestamp>` first. Text that was edited by hand is left alone; the summary lists any text box that still describes separate dot plots. A second run finds nothing to merge.

## Layout and thumbnails

The question views and `Cell atlas` use fixed pixel positions with about 25 px between charts, so the page can be scrolled from the gaps without the mouse wheel zooming a chart. Each view keeps its width; it gets a little taller.

The scripts don't write thumbnails (`viewImage`). Open each new view with edit permission and save it (Ctrl+S) before copying `views.json` to a server where users only have view permission. A viewer who opens a project whose first view has no thumbnail gets a "Forbidden" save error.
