# AnnData → MDV mapping

How a Scanpy/AnnData object (`.h5ad`) is represented once it has been imported
into an MDV project. This is the conceptual bridge between the two data models;
for the on-the-wire column format see [datasource](./datasource.md) and for the
expression-matrix link see [observable column links](<./observable column links.md>).

## Two models, one projection

AnnData is a container of **typed slots** — `X`, `obs`, `var`, `obsm`, `varm`,
`uns` — each with its own semantics, all carried by a single in-memory object and
frozen to disk by `write_h5ad`. It is optimised for *analysis*: you keep mutating
one object and it accumulates the history of everything you've done.

MDV is the opposite. It has exactly **one structural primitive** — the
*datasource* (a flat table of index-aligned columns, backed by one group in
`datafile.h5` and one entry in `datasources.json`) — plus **one relational
primitive** — the `rows_as_columns` *link* that joins two datasources as a matrix.
It is optimised for *interactive visualisation at scale*: every column is an
independently loadable, indexed, quantile-summarised dataset that a chart can bind
to.

Importing an AnnData object therefore **projects** its typed slots onto those two
primitives. Nothing in MDV is "a UMAP" or "a DGE result" — everything is a column,
a matrix link, or a whole datasource.

## The slot mapping

| AnnData slot | What it is | MDV home |
|---|---|---|
| `adata.obs` | per-cell metadata (cluster labels, QC) | columns on the `cells` datasource |
| `adata.var` | per-gene metadata | columns on the `genes` datasource |
| `adata.X` / `adata.layers[...]` | expression matrix (cells × genes) | the `links.rows_as_columns` matrix bridge (`cells/gene_scores`, the `gs` subgroup) |
| `adata.obsm` | dimensionality reductions (cells × k) | **flattened** into scalar columns on `cells` (`X_umap_1`, `X_umap_2`, `X_pca_1..3`, ...) |
| `adata.varm` | per-gene loadings (genes × k) | flattened into columns on `genes` (`PCs_1..3`) |
| `adata.uns` | unstructured grab-bag (incl. DGE results) | **— no slot —** see below |
| `adata.raw` | pre-HVG full expression | **dropped on import** (see caveats) |

## How a UMAP is represented

There is no UMAP type in MDV. On import, `_add_dims` (`conversions.py`) takes each
`obsm` entry of shape `cells × k` and unrolls it into `k` ordinary `double`
columns named `{name}_{i+1}` — so `obsm['X_umap']` becomes two columns,
`X_umap_1` and `X_umap_2`, type-identical to `n_counts` or `percent_mito`.

The "UMAP-ness" then lives entirely in:

1. **the chart config** — a scatterplot whose x = `X_umap_1`, y = `X_umap_2`,
   coloured by e.g. `leiden`; and
2. **the naming convention** — the `X_umap_*` prefix.

The data model does not know those two columns are an embedding rather than two
arbitrary metrics. The same is true of PCA (`X_pca_*`), t-SNE (`X_tsne_*`), etc.
`max_dims` (default 3) caps how many components of each reduction are imported.

## The `uns` gap — where job outputs go

`adata.uns` is the one slot with no MDV equivalent. It is AnnData's unstructured
dict: DGE results, colour palettes, run parameters — anything that doesn't fit the
`obs`/`var`/`X`/`obsm` axes. Scanpy writes differential-expression results there
(`sc.tl.rank_genes_groups` → `uns['rank_genes_groups']`), **not** into `var` —
because a DGE statistic is not a function of the gene alone; it is keyed by the
contrast (which two cell groups), the cell selection, and the test method.

For a *visualisation* tool an untyped dict is useless — you cannot bind a chart to
a nested record array. So MDV deliberately has no `uns`. The faithful translation
of "a `uns` entry you want to visualise and slice" is to **promote it onto MDV's two
primitives** — a column on an existing datasource, a `rows_as_columns` matrix, or a
whole new datasource — chosen by the output's **shape**. That shape-routing is the job
framework's output contract; see
[ADR-0006](../adr/0006-one-data-only-tool-spec-job-keyed-idempotent-output.md). The
POC's first tool, `concat_columns`, takes the simplest route: one new column on the
chosen datasource.

> **DGE is deferred.** A DGE result is keyed by *(gene, comparison)* — the rich
> `matrix` + `new_entity` case (a `contrasts` datasource + `genes × contrasts`
> matrices, able to hold genes outside `var` on its own index). That design is parked
> with DGE: see [ADR-0003](../adr/0003-read-expression-from-mdv-store-not-source-h5ad.md)
> (deferred).

## Where this lives on disk

An MDV project directory holds the manifest and the data side by side:

```
project/
  datasources.json   # the manifest: array of { name, size, columns[], links?, images? }
  datafile.h5        # the data: one HDF5 group per datasource
  views.json         # saved chart layouts
  state.json         # project state / permissions
```

`datasources.json` carries **no values** — it is a schema/index. Each column's
`field` is the key into the datasource's h5 group:

- a scalar column `cells/n_genes` is a 1-D array of length `size`;
- a `text16` column is stored as `uint16` *codes*; the `values` array in the
  manifest is the decode table;
- the expression matrix is a subgroup `cells/gene_scores` holding a flattened `x`
  array (`n_cells × n_genes`) plus a `length` scalar — exposed through the
  `rows_as_columns` link, not as a column.

A new datasource (e.g. a job output) is born via `add_datasource`
(`mdvproject.py`), which creates the h5 group, writes the columns
(`add_column_to_group`), and appends the manifest entry
(`set_datasource_metadata`). It is removed cleanly via `delete_datasource`
(`del h5[name]` + drop the manifest entry + drop dependent views).

## Frontend load

`datasources.json` is fetched once at bootstrap (`src/modules/projectRuntime.ts`)
and `ChartManager` (`src/charts/ChartManager.js`) instantiates one `DataStore` per
manifest entry. A datasource added *after* bootstrap is therefore not visible until
the next load — which is why a freshly written job output currently relies on the
existing `needs_refresh` → reload path.

## Code references

| Concern | Location |
|---|---|
| AnnData → MDV import (entry point) | `python/mdvtools/conversions.py` · `convert_scanpy_to_mdv()` |
| `obsm`/`varm` flattening | `conversions.py` · `_add_dims()` |
| `X` → CSC store / compute-at-import | `conversions.py` · `_coerce_x_to_csc_without_nonfinite()`, `_prepare_x_umap_and_leiden()` |
| read/write the manifest | `python/mdvtools/mdvproject.py` · `datasources` property |
| create / delete a datasource | `mdvproject.py` · `add_datasource()`, `delete_datasource()` |
| write a column + deduce its metadata | `mdvproject.py` · `add_column_to_group()`, `get_column_info()` |
| column / datatype format spec | [datasource](./datasource.md) |
| the expression-matrix link | [observable column links](<./observable column links.md>) |

## Caveats

- **Slot semantics are not preserved.** Once `obsm['X_umap']` becomes
  `X_umap_1/2`, MDV only knows they are an embedding by naming convention and chart
  usage. A clean MDV → AnnData round-trip has to *re-infer* slot membership.
- **`adata.raw` is dropped on import.** DGE is meant to run against the full
  expression set, but the `genes` datasource is built from `var` (the HVG subset),
  so its gene set is usually smaller than `.raw`. Whether to preserve `.raw` is OQ1,
  **deferred with DGE** — it does not affect the `concat_columns` POC.
- **`uns` does not round-trip** — there is no slot to put it back into; it is
  reconstructed from whatever datasources/columns it was promoted into.

## External references

- AnnData data model — https://anndata.readthedocs.io/en/latest/
- `.h5ad` / Zarr on-disk layout — https://anndata.readthedocs.io/en/latest/fileformat-prose.html
- Scanpy `rank_genes_groups` (DGE → `uns`) — https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html
