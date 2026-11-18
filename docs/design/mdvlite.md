# `mdvlite`: serve MDV over SpatialData with no project and no data stack

> **Status:** design note, not a decision. It follows from
> [ADR-0023](../adr/0023-write-project-data-to-icechunk-and-retire-h5.md) and leans on
> [ADR-0024](../adr/0024-one-column-type-contract-built-on-zarr-data-types.md).

## The idea

```bash
uvx mdvlite serve ~/data/xenium.zarr
uvx mdvlite serve https://example.org/stores/visium.zarr
```

A small Python package that opens MDV's charts, with linked filtering, over SpatialData stores,
local or remote. It needs no existing MDV project and no conversion step. It does for MDV what
[`sdjs`](https://github.com/Taylor-CCB-Group/SpatialData.js/tree/main/python/sdjs) does for the
SpatialData.js viewer: `sdjs` has no dependencies at all (a standard-library server plus the viewer
bundled as static files), finds `*.zarr` stores under a path, and serves them with a picker.

What the user builds is saved in a form a full MDV deployment can later take in as a project.

## The earlier lite build, and what to do differently

`mdvtools_lite_build/` (#172) published a lighter `mdvtools` by cherry-picking a few modules with
Hatchling `force-include`. Its dependencies were not locked, and the include list was brittle:
the published wheel was not what an editable install ran. ADR-0001 removed it. `mdvlite` should
avoid all three problems:

- **Its own source tree, with nothing copied or imported from `mdvtools`.** That is ADR-0001's
  "pattern B" (genuinely separate packages), not the subset build ADR-0001 rejected.
- **Locked dependencies** (`uv.lock`), and as few of them as possible. Ideally none, like `sdjs`.
- **Anything both packages need is shared as a contract, not as code.** The column-type table
  (`column_types.json`, ADR-0024) and the project descriptor below are data both sides read.

## Why this only becomes neat once h5 goes

Today "an MDV project" means `datafile.h5` plus a `datasources.json` describing it, produced by
`convert-spatial`, which needs the whole Python data stack. `mdvtools`' slim core still installs
around 1.2 GB of dependencies (ADR-0001), and `serverlite.py` is built on `MDVProject` and h5.

Under ADR-0023 the source store *is* the data, and the browser reads columns from it through a
zarr `DataLoader` (the deployment-swappable seam the front end was designed around). The server
is left to serve the front end and the stores' files, and to persist a few small JSON files.
None of that needs h5, pandas or a codec library.

## A project descriptor instead of `datasources.json`

Rather than persisting `datasources.json` in its current shape, `mdvlite` persists a small
**project descriptor**: the list of input stores and a formal description of how they are
interpreted. It amounts to the arguments `convert-spatial` would have been given. The browser
manifests the `datasources.json`-shaped config the front end expects from the descriptor plus
the stores' own metadata. The descriptor is the source of truth; the manifested config is derived
and can be cached.

`convert-spatial`'s options fall into four groups, and only two belong in the descriptor:

| Group | Options | In the descriptor? |
|---|---|---|
| **Which stores, and where they are** | the store paths, `--link`, `--spatial-url`, `--batch` | Yes: each input store, and whether it is copied, linked or remote |
| **How they are interpreted** | `--table-handling`, `--obs-datasource-name`, `--var-datasource-name`, `--link-name-column`, `--point-transform` | Yes: this is what makes the result reproducible |
| **Rewriting the store for the browser** | `--browser-preset`, `--pyramids`, `--image-codec`, `--image-preset`, `--csc-tables`, `--morton-points`, `--points-zarr` | No: an offline step that produces a source store. The descriptor may record that it was applied |
| **Computing new data** | `--compute-x-umap`, `--leiden-resolution` | No: these produce new columns, which are Project columns (ADR-0023) or the output of a job |

`--output_geojson` belongs in none of them: it writes region GeoJSON that only legacy charts read,
and is being retired (#670). The rest (`--density`, `--serve`, `--preserve-existing`,
`--verbose`) still need sorting into these groups.

**The equivalence contract.** For the same stores and the same descriptor, the manifested
datasources should match what `convert-spatial` writes, apart from fields that only exist
because of h5. This is testable: run both over the synthetic SpatialData fixtures and compare.
Once it holds, `convert-spatial` itself could write the descriptor, and full projects and
`mdvlite` would share one format. `datasources.json` would then be derived for both.

**`state.json` too.** It holds the view list and initial view (`all_views`, `initial_view`),
which belong to the project, alongside `permission` and the display name, which belong to the
deployment serving it. The display name is there as the catalog's recovery copy. `mdvlite`
persists the view list and initial view with `views.json`. On import, a full deployment keeps
those and sets the permission and name itself.

## Serving, reads and writes

- **Serving.** Reuse what `sdjs` already does: store discovery, `catalog.json`, tokenised URLs,
  range requests. Either depend on `sdjs` or share its server module.
- **Reads.** All column data is read in the browser from the source stores. A remote store is read
  from its own URL, which needs CORS on the remote. There is no `/get_data`.
- **Datasources.** One per table (ADR-0023 decision 4), named deterministically from store and
  table, because views bind to datasources by name and an import must reproduce the names.
- **Writes.** The descriptor, `views.json` and the project part of `state.json`, to a sidecar
  directory next to the stores or one the user names. Column edits, annotations and jobs are out
  of scope at first. With the `icechunk` package as an optional extra, a sidecar Project
  repository could add them later.
- **Import into full MDV.** A deployment creates a project from the descriptor (copying or linking
  the stores as it says) and takes `views.json` and the view list.

## The front-end bundle

MDV's built front end is about 59 MB. A package with no (or very few) Python dependencies and a
59 MB bundle is still light next to today's install, so this is not a blocker. Later, the bundle
could be split, and a more formal monorepo layout could separate catalog, spatial, genome
browsing and so on. Trimming legacy charts would also help.

## Open questions

- **The descriptor's name and schema,** and whether it lives beside `views.json` or replaces
  `datasources.json` in full projects too.
- **Where the manifest is computed.** The browser is the natural place: it already reads category
  values and statistics through zarrita and SpatialData.js, and doing it there keeps Python
  dependency-free. It is also what full MDV's runtime-derived datasources would use
  (spatial-tables theme 4). Does the server cache the result?
- **Non-SpatialData zarr** (ADR-0023 decision 1): AnnData zarr, and plain zarr groups of 1-D
  arrays. Which layouts does the descriptor describe first?
- **Several tables per store.** It's one datasource per table, but what does the first view show,
  and how do linked tables filter across each other before element relations (#555) land?
