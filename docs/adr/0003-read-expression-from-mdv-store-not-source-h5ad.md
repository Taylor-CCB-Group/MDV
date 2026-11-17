# Read the expression matrix from MDV's CSC store, not the source `.h5ad`

**Status:** Deferred — DGE-specific; out of scope for the first jobs-framework POC. The POC's
first tool is a trivial concat-columns job, which needs none of this (no expression matrix, no
CSC read). Revisit when DGE work resumes. Prior art: PR #352.

**Amended** by ADR-0023 (proposed, Taylor-CCB-Group/MDV#665): see *Amendment: Icechunk project storage* at the end.

> Parked 2026-06-11: the jobs-framework POC pivoted away from DGE-first to build the framework
> against a throwaway concat-columns tool. This decision — along with ADR-0006's DGE worked
> example, the Stage-3 `contrasts`/matrix design, and OQ1 — is deferred with DGE, not discarded.

## Context & decision

DGE needs a genes×cells matrix as input. Two sources are possible: re-read the original
uploaded `.h5ad`, or read MDV's own store. We **read from the store** — the rows-as-columns
subgroups under the `cells` datasource (`gs` holds `adata.X`; each named layer is its own
subgroup; see `conversions.py` `add_rows_as_columns_subgroup`).

The slice the worker runs on is the **cells in the current client-resolved selection**
(the resolved index list the frontend sends) × **all genes**. It is **not** a re-read of the
`.h5ad`: the selection indexes the `cells` datasource rows, which are 1:1 with the stored
matrix rows.

## Why

- **No drift** — DGE runs on exactly the values in MDV's store, including edited `obs` /
  renamed labels.
- **No cell-identity join** — matrix rows *are* the `cells` datasource rows, in the same
  order; there is no fragile mapping back to an external h5ad.
- **Labels from the store, never the client** — for the client's resolved selection (*which*
  cells), the group **labels are re-read from the store**, not trusted from the client.
  Selecting cells is the client's job; *what those cells are* (their `obs` labels) is the
  store's. (A corollary of no-drift, and the load-bearing half of the D3 selection contract.)

## The bulk reader (required consequence)

A new **bulk** CSC reader is needed. The store is genes-as-columns: `p` is indexed over
genes, `i` holds cell-row indices (`mdvproject.py` `_read_subgroup_matrix_column`). The bulk
path reads `p`/`i`/`x` **once** → `scipy.sparse.csc_matrix((x, i, p), shape=(n_cells,
n_genes))` → AnnData.

Do **not** reuse the per-gene path (`get_datasource_as_dataframe` →
`_read_wrapper_expression_column`): it opens/closes the h5 per gene and returns Python
`list[float]` — O(n_genes) and pathologically slow for a full matrix.

## Consequences

- A per-cell subset is CSC's **expensive** axis. Row-slicing (selecting cells) on a CSC is
  slow, so the reader materializes the **full** matrix (all `nnz`) and subsets cells in
  memory (convert to CSR first). Gene-subsetting is the cheap axis and DGE doesn't use it (it
  needs every gene). Memory ≈ O(full nnz): fine at PBMC3k scale, and this is **the** scaling
  pressure point — the same reason the per-gene path was rejected, now on the bulk path.
- "Exactly what the user sees" is scoped to **cell membership + labels + `obs` edits** — the
  things visible and editable in the UI. It does **not** cover the matrix *representation*
  (raw counts vs log-normalized), which never appears on screen. Representation-correctness
  is **OQ1** (open — confirming whether imported `adata.X` is raw or log-normalized), not a
  guarantee of this decision. `rank_genes_groups(method="wilcoxon")` conventionally expects
  log-normalized input, so OQ1 is a correctness boundary, handled there (an expression-layer
  choice recorded in provenance — the validation guard was **dropped**; we trust the input),
  not here.

## Amendment: Icechunk project storage (ADR-0023, proposed)

ADR-0023 ([#665](https://github.com/Taylor-CCB-Group/MDV/pull/665)) retires `datafile.h5`. Bulk data
stays in read-only **source stores**, and everything MDV writes goes to the project's **Icechunk
repository**. When DGE resumes, this decision holds and its reader changes:

- **"The store" is wherever the matrix lives.** For data imported from AnnData or SpatialData that is
  `X` (or a layer) in the source store. For a project upgraded from h5 it is the repository. Either
  way it is still MDV's data and never a separate `.h5ad`, so the no-drift argument is unchanged.
- **The bulk reader reads the sparse layout the source wrote**, instead of MDV's gene-column CSC
  subgroups. If the source is CSR (cells as rows, as AnnData usually writes it), selecting cells
  becomes the cheap axis, and the full-`nnz` pressure point in *Consequences* goes away for those
  stores.
- **Labels are resolved the way `/get_data` resolves them.** An `obs` column may come from the source
  store or be a project column holding the user's edits, and edits win. The materializer has to use
  the same column routing as the rest of the server, or the no-drift guarantee breaks for edited
  labels.
