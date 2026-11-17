# Write project data to an Icechunk repository and retire `datafile.h5`

**Status:** proposed. The jobs ADRs on `feature/dge_poc_local` (0004, 0012, 0014) assume HDF5 and
are to be amended alongside this before that branch merges. Numbered 0023 because branches
already use 0004–0022.

## Context

Every column MDV shows is stored in `datafile.h5`. That includes SpatialData tables, which
`convert-spatial` also leaves in the project's copy of the zarr store. The spatial viewer reads
the columns it finds in the zarr table straight from zarr, and the DataStore reads everything else
from h5 over `/get_data`; the two are joined by key
([table_association.ts](../../src/react/spatialdata/table_association.ts)). The same table
exists twice, with no rule for which copy wins.

h5 costs us in other places too:

- **Locking.** HDF5 refuses to open a file for writing in a process that already has it open for
  reading. The project's fcntl lock does not exclude threads, hot paths do not take it (#538), and
  conversions take no lock at all (#623). ADR-0014 exists largely to route around this.
- **Merged datasources.** When several similar tables become one datasource, every column is
  materialised for the union of their rows, so loading and filtering a column pays for rows it
  has nothing to do with.
- **History.** The project JSON history utility (#557, #563) is a versioning layer we built by
  hand.
- **Types.** Values are narrowed on the way in (ADR-0024).

The notes in `docs/design/spatial-tables/` explored reading tables from zarr in the browser while
keeping every write in h5. This ADR replaces that split.

## Decision

1. **Source stores are read and never written.** Bulk data stays in zarr stores MDV does not
   modify: SpatialData stores, AnnData zarr, and zarr layouts that are not SpatialData at all. The
   server keeps serving them as files through `/spatial/<path>`, and the browser reads them
   directly, as it does today.
2. **Everything MDV writes goes to the project's Icechunk repository.** There is one repository
   per project, so one commit can cover columns on several tables. It holds:
   - project columns: user edits and annotation columns, job outputs, and columns imported
     without a source store;
   - the history of the project JSON (`datasources.json`, `views.json`, `state.json`), replacing
     the snapshot utility from #557;
   - later, auxiliary spatial stores: shape and label annotations, and editable scenes that carry
     user-written coordinate-system transforms. These are where mutable SpatialData gets
     prototyped, ahead of the upstream work it needs.
3. **The server owns the repository.** For now browsers neither read nor write Icechunk. The
   server reads project columns from the repository and serves them through `/get_data`, and
   every write goes through a server endpoint, which is where permissions are checked. Keeping
   reads there too means one reader to build first; browser reads are a later option (see
   *Considered options*). The server serialises commits to a project's repository, because
   Icechunk documents its local-filesystem storage as unsafe under concurrent commits: both can
   report success while one is lost.
   The project write lock is an fcntl lock, which excludes other processes but not threads or
   greenlets in the same process, so commits made inside the server process also hold an
   in-process lock.
4. **A datasource is one table.** Its columns load from wherever they live, a source store or the
   repository, and are not copied into a merged datasource. Similar tables are related through
   element relations (#555) and filtered across, not concatenated. `convert-spatial`'s by-region
   merging carries on for now and is not built on further.
5. **h5 is retired for active projects.** An upgrade script turns an h5 project into a repository;
   every column moves into it, since there is no source store behind it. A downgrade script goes
   back to h5 so a project can be opened by an older MDV: it writes source columns into h5,
   narrows types the way h5 projects always have, and drops what h5 cannot hold (history,
   auxiliary stores). The server recognises which kind of project it has opened and keeps serving
   h5 projects through the transition.
6. **Static export writes plain zarr.** It flattens the repository's current snapshot, with the
   source stores, into zarr that the static build reads directly, replacing the `.gz` files and
   their `.json` byte-range index. The static build needs no Icechunk reader.

ChatMDV's use of `.h5ad` and CSV files is out of scope and handled separately.

## Considered options

- **Keep h5 behind a storage interface.** Wraps the thing we want to remove and keeps its locking
  and narrowing.
- **Read zarr in the browser, keep writing h5** (the split in the design notes). Leaves two copies
  of each table and keeps every h5 problem on the write side.
- **Plain zarr for writes.** No transactions or history, so concurrent writes and JSON history
  would need machinery of our own.
- **Write into the source SpatialData store.** Source stores can be shared, linked (`--link`) or
  remote (`spatial.url`), and SpatialData has no mutation model yet.
- **Reading Icechunk in the browser.** Deferred rather than rejected. A v2 repository can be read
  from a plain static file server without listing, and the community `icechunk-js` package
  (read-only, built on zarrita) opens one from a URL, so MDV could later serve a repository the
  way it serves source stores. The official `@earthmover/icechunk` package documents browser
  access only for S3-compatible hosting. Writes stay on the server either way.

## Consequences

- Column metadata records where a column's values live (a source store and path, or the
  repository), and the DataLoader routes on it. This is the `backing` hint from the design notes,
  generalised beyond h5 versus zarr.
- `icechunk` becomes a Python dependency. zarr is already on v3.
- **Where a repository lives.** By default, local-filesystem storage in the project folder, which
  needs no S3 and works on macOS, Linux, and Windows through Docker or WSL2 (MDV's documented
  install). Icechunk's local-filesystem storage does not work on native Windows, so a native
  Windows install of `mdvtools` cannot open a repository; the downgrade script is its fallback.
- **Serialising commits depends on the lock reaching every committer.** That holds for one server
  on one machine. On network filesystems (NFS, Lustre, GPFS, SMB), which Icechunk's docs do not
  cover, both Icechunk and the fcntl lock need testing before a project lives there. Job workers
  never commit (ADR-0004), so on a cluster only the server's node writes.
- **Shared projects need object storage.** Two independent deployments sharing a Project root by
  mount or bucket sync could commit to the same repository through a filesystem, and one commit
  could be lost. Such deployments put the repository on object storage with conditional writes
  (S3, GCS, Azure, MinIO) instead.
- The fcntl write lock stays, now serialising commits, alongside an in-process lock. HDF5's
  same-process restriction goes. The jobs ADRs on `feature/dge_poc_local` carry amendments for
  this: ADR-0014's gevent reason for running owner steps in a subprocess still holds and its HDF5
  reason does not; ADR-0004's trays can stay HDF5, since they are scratch, but materialize reads a
  snapshot pinned at submit rather than taking the read lock.
- Two project formats coexist until upgrades are done, so tests need fixtures of both.
- Per-table datasources lean on cross-datasource views and filtering (design themes 4 and 5).
