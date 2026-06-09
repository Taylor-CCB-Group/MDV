# Read the expression matrix from MDV's CSC store, not the source `.h5ad`

**Status:** accepted — POC design decision (DGE jobs framework); not yet implemented on
this branch (prior art: PR #352).

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
  log-normalized input, so OQ1 is a correctness boundary, handled there (layer dropdown +
  provenance + a raw-counts guard), not here.
