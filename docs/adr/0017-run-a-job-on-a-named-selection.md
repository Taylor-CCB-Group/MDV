# Run a job on a named selection of rows

**Status:** accepted (jobs framework).

## Context

A user who has a UMAP embedding plots `UMAP_1` against `UMAP_2`, lasso-selects a group of cells and
wants to run UMAP again on only those cells. They also want to name the group, such as "T cells", so
they can tell the subset and its results apart later. Until now every job ran on all rows of its
datasource.

MDV filters a datasource by combining every active filter: a lasso, range and category filters, and
gates. The bar shows how many rows pass (`DataStore.filterSize`), and `getFilteredIndices()` returns
them. Filters change freely after a job is submitted, so the job has to fix the exact rows it was
given, and restart recovery (ADR-0005) has to see the same rows.

## Decision

A job can run on the datasource's current filtered rows, whatever filters produced them. The user
names the selection, the owner fixes the rows at submit, the materializer cuts them out of the input,
and ingest writes the result back into full-length columns with NaN in the rows outside the
selection.

### How the rows reach the server

The browser sends the sorted row indices in a top-level `rows` field of `POST /jobs`, next to
`params`. At submit the manager checks them and writes them to a sidecar file,
`<project>/jobs/selections/<job_id>.npy`, before it writes the job record. The record holds
`input_filter_hash`, a hash of the indices, and `selection_rows`, their count. It never holds the
indices themselves.

We keep the indices out of `params` because the record goes into every `GET /jobs` response and is
piped to each owner step on stdin (ADR-0014), so a list of a million indices would travel with every
poll of the job list. We send indices rather than writing the selection to a column first, because a
column made in the browser reaches the server only when the user saves the view, and a column per run
would clutter the project.

The sidecar is part of the job's durable history and lives as long as its record. It costs 4 bytes
per selected row, so 4 MB for a million cells.

The manager rejects `rows` that are empty, are not whole numbers, fall outside the datasource's size
or contain duplicates. It also rejects `rows` sent to a tool that cannot run on a selection, and a
selection name sent without rows or rows sent without a name. A selection that is too small for the
tool, such as fewer cells than `n_neighbors`, is not the gate's concern (ADR-0006): the worker fails
the job with its own error.

### How a tool declares it

`ParamSpec` gains a `selection` type, bound to its datasource through `options_from`. Its value in
`params` is the selection's name, and it is empty when the job runs on all rows. A tool that lists a
`selection` param can run on a selection, and the form renders it as the choice between all rows and
the current selection, together with the name.

Only the matrix input shape supports a selection for now, so UMAP is the only tool that declares one.
The columns input shape gets it when a tool needs it, which is also when we will choose a missing
value for text columns.

### Where the subset is applied

The owner's materializer reads the sidecar and writes only the selected rows into the tray, as
ADR-0004 already described for DGE's "selected cells × all genes" tray. Workers do not change and do
not know a subset exists. For a sparse matrix the materializer keeps the CSC entries whose row is
selected and renumbers the rows. For a dense matrix it takes the selected rows from each gene's slice.
The manager passes the sidecar's path to the materialize and ingest owner steps as an optional fourth
argument, so the steps still need no `records_root` (ADR-0014).

### What the output looks like

The selection's name becomes the default `output_name`: "T cells" gives `UMAP_T_cells`, so the
columns are `UMAP_T_cells_1` and `UMAP_T_cells_2`. The user can still edit `output_name`. Running
again with the same name replaces those columns in place, and the newest job owns them (ADR-0009).

Ingest reads the sidecar, places the worker's result in the selected rows and writes NaN in every
other row, so output columns stay full length. MDV already treats NaN as missing in a double column:
the DataStore shows it as "missing", and `minMax` and the quantiles skip it.

### Provenance

The provenance pointer on each output column gains `selection: {name, n_rows}` when the job ran on a
selection, so a column can describe itself as "UMAP on T cells (2,314 cells)" without opening the
record. `content_hash` already includes `input_filter_hash`, so two runs on different rows count as
different analyses.

### The dialog and the job list

When a filter is active on the chosen datasource, the Run analysis dialog defaults to "Current
selection (2,314 of 2,700 cells)" and offers all cells as the alternative. Submit stays disabled
until the selection has a name. With no filter active the dialog looks as it did before. The job list
gains a "Cells" column that shows the selection's name and count, or "All".

## Considered options

- **Indices inside `params`.** Simple, but they would travel with every job list poll and every
  owner step.
- **Write the selection to a column and send its name.** Needs a forced view save or a new endpoint,
  and leaves a column in the project for every run.
- **Give the worker an index list.** Every worker would have to handle subsets, and the tray would
  stay full size.
- **A category column marking the selected rows.** Gates already do this, so the name stays a label
  and part of the output names.
- **Default to all rows when a filter is active.** Safer against a forgotten filter, but it adds a
  step to the workflow this feature is for. The count shown in the dialog and the different default
  output name guard against a forgotten filter instead.

## Consequences

- Records written before this change have no `selection_rows` and no sidecar, and run on all rows.
- The WGL scatter plot has no NaN guard, so a NaN coordinate reaches the shader. GPUs normally drop
  such a point, but GLSL leaves the behaviour undefined, so this has to be checked in the browser.
- A WGL scatter plot does not redraw when a re-run replaces its columns in place, because it keeps
  the typed arrays it was built with. Re-running a named selection hits this known issue.
- Row indices assume the datasource's rows never change order or count after the job is submitted,
  which holds for MDV datasources today.
