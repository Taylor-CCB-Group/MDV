# `mdvlite`: serve MDV over SpatialData with no project and no data stack

> **Status:** design note, not a decision. It follows from
> [ADR-0023](../adr/0023-write-project-data-to-icechunk-and-retire-h5.md) and leans on
> [ADR-0024](../adr/0024-one-column-type-contract-built-on-zarr-data-types.md).

## The idea

```bash
uvx mdvlite serve ~/data/xenium.zarr
uvx mdvlite serve https://example.org/stores/visium.zarr
```

A small Python package that opens MDV's charts, with linked filtering, over a SpatialData store,
local or remote. It needs no MDV project and no conversion step. It does for MDV what
[`sdjs`](https://github.com/Taylor-CCB-Group/SpatialData.js/tree/main/python/sdjs) does for the
SpatialData.js viewer: `sdjs` has no dependencies at all (a standard-library server plus the viewer
bundled as static files), finds `*.zarr` stores under a path, and serves them with a picker.

Views the user builds are saved, as `views.json`, in a form a full MDV deployment can later
import as a project.

## Why this only becomes neat once h5 goes

Today "an MDV project" means `datafile.h5` plus `datasources.json` describing it, produced by
`convert-spatial`, which needs the whole Python data stack (spatialdata, scanpy, pandas, h5py).
`mdvtools`' slim core still installs around 1.2 GB of dependencies (ADR-0001), and
`serverlite.py` is built on `MDVProject` and h5. A lite server on that basis would either drag
the stack along or need its own converter.

Under ADR-0023 the source store *is* the data. The browser reads columns from it through a
zarr `DataLoader`, which is the deployment-swappable seam the front end was designed around.
What is left for the server is small:

- serve the bundled front end, and the store's files with range requests;
- provide a `datasources.json` describing the store's tables;
- persist `views.json` (and `state.json`).

None of that needs h5, pandas or a codec library.

## Sketch

- **Package.** Its own PyPI package, `mdvlite`, with its own source tree and no dependency on
  `mdvtools`. That is ADR-0001's "pattern B" (genuinely separate packages), not the subset build
  ADR-0001 removed. Dependencies: none, or as few as `sdjs` has.
- **Serving.** Reuse what `sdjs` already does: store discovery, `catalog.json`, tokenised URLs,
  range requests. Either depend on `sdjs` or share its server module.
- **Reads.** All column data is read in the browser from the source store. A remote store is read
  from its own URL, which needs CORS on the remote. There is no `/get_data`.
- **Datasources.** One per table (ADR-0023 decision 4), named deterministically from store and
  table, because views bind to datasources by name. Each column's metadata comes from the
  column-type contract (ADR-0024).
- **Writes.** Views and state only, to a sidecar directory next to the store or one the user
  names. Column edits, annotations and jobs are out of scope at first. With the `icechunk`
  package as an optional extra, a sidecar Project repository (ADR-0023) could add them later.
- **Import into full MDV.** A deployment creates a project whose source store is the same store,
  copied or linked, and takes the sidecar `views.json`. That works only if datasource names come
  out the same in both, which is one more reason they must be deterministic.

## Open questions

- **Where is `datasources.json` manifested?** In Python, reading the store's zarr JSON metadata is
  standard library, but category values and column statistics live in (usually compressed)
  arrays, so Python would need codecs. The browser already reads them through zarrita and
  SpatialData.js. Manifesting in the browser keeps the Python side dependency-free, and it is the
  same code full MDV's runtime-derived datasources would use (spatial-tables theme 4). The server
  could then cache the result.
- **Bundle size.** MDV's built front end is about 59 MB, against `sdjs`'s much smaller viewer.
  Is a slimmer build with only the charts `mdvlite` needs worth having?
- **Non-SpatialData zarr** (ADR-0023 decision 1): AnnData zarr, and plain zarr groups of 1-D
  arrays. Which layouts does `mdvlite` recognise first?
- **Several tables per store.** It's one datasource per table, but what does the first view show,
  and how do linked tables filter across each other before element relations (#555) land?
