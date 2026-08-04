# Remote executors assume a shared filesystem at a matching path; data movement is a deferred seam, not part of the executor

**Status:** accepted — POC design decision (jobs framework); remote-staging and pull-agent backends deferred.

## Context & decision

The next backend is a `SlurmExecutor` (ADR-0008 blesses it as a peer class behind
`submit`/`poll`/`locate_result`). Writing it surfaced a question the POC never had to answer
while everything ran as a local subprocess: **how does the tray reach the worker, and the output
come back, when the worker runs somewhere else?**

The decision: **the first `SlurmExecutor` — and any near-term remote executor — assumes a shared
POSIX filesystem visible at a *matching path* on both the owner and the compute node.** Under that
assumption the courier model (ADR-0004) is unchanged: the owner materializes the tray into the
workspace, the worker reads and writes that same directory by path, the STATUS marker is visible to
the owner the instant it appears, and `locate_result` is `workspace/output` exactly as it is
locally. The executor stays the thin CLI class ADR-0008 describes; it moves **no bytes**.

This ADR records that precondition, why lifting it (a *truly* remote owner with no shared FS) is a
separate subsystem rather than a config value, and the evolution ladder for that case. It builds on
ADR-0004 (which already names staging as a *separate seam* and cites Pulsar/Bazel) and ADR-0007
(workspace = ephemeral scratch on `$SCRATCH`/shared FS); it does not restate them. It is distinct
from the coming resource/placement ADR (`--mem`/`--cpus`/`--gres`), which is about *which node* a
job asks for, not *how data reaches it*.

## The precondition: a shared FS at a matching path

"Shared filesystem" is the real discriminator, not "local vs remote" (ADR-0004). A remote owner is
fine **as long as it mounts the cluster's filesystem at the same path the compute nodes use.** That
is achievable off a single machine — NFS / Lustre / GPFS / BeeGFS, or a managed cloud FS (FSx for
Lustre, EFS, Filestore) — and is the ordinary HPC/analysis-server setup. When you can arrange it,
the whole "remote" problem dissolves back into this design.

Two things it demands, both easy to miss:

- **Login/SSH access is a weaker grant than a mountable export.** "The owner has read/write access
  to the cluster" usually means a login shell — *not* that cluster storage is exported to the
  owner's host as a mountable POSIX filesystem. The second is an admin + networking decision
  (NFSv4/Kerberos ACLs, low-latency line of sight; NFS over a WAN is painful). Arrange it; don't
  assume it.
- **Paths must line up.** The workspace's absolute path is carried in the worker argv and Slurm
  `--chdir`. If the owner sees `/mnt/cluster/scratch/<job_id>` while the node sees
  `/scratch/<job_id>`, the handoff breaks even though the bytes are shared. Mount at the **identical
  path** on both sides (or add a path-translation shim — more machinery). Same-path mounting keeps
  the executor trivial.

## Why lifting it is a subsystem, not a config value

If the owner genuinely cannot mount the cluster FS, the framework must copy the tray onto the
cluster before `sbatch` and copy the output back after — "network staging." This is done in
practice (Pulsar, Nextflow's file-system abstraction, Globus), but it is not a flag:

- **It re-adds the coupling the courier model deleted.** Path-on-shared-FS is zero-copy and needs no
  protocol. Staging needs a transfer mechanism, credentials, integrity checks, retry/resume, and a
  way to *know the upload finished* before submit. That is a new service, not a setting.
- **It hits the largest data.** The jobs pushed to Slurm are heavy *because* the tray is a multi-GB
  sparse matrix. The shared-FS materializer does a sparse copy on the same disk and never moves
  those bytes; staging ships GBs per job, twice, competing with the web process's bandwidth. Even
  on-cluster, HPC guidance already stages home/project → fast per-job `/scratch` to dodge NFS
  congestion — so staging is a real tier, which is exactly why it deserves to be its own seam.
- **It inverts the completion signal.** Today the STATUS marker is primary because the owner shares
  the FS; `poll` is secondary (ADR-0008). With no shared FS the owner cannot cheaply `stat` the
  marker, so `poll` (`sacct`/REST) becomes primary and the marker is staged. That is a signal-lead
  flip, anticipated by the existing primary/secondary split — not a redesign.
- **It muddies durability.** ADR-0005 recovery assumes the workspace is on disk where the owner
  finds it. If the real workspace lives on the cluster and the owner holds copies, "source of truth
  after a crash mid-transfer" becomes a live question.

## Evolution ladder for the no-shared-FS case

Ordered by how much they disturb this design. The sealed worker (ADR-0004) is what keeps the first
three *additive* — because the worker only touches a local workspace path, a data-movement layer
slides underneath it without the worker knowing, exactly as it never knows it is on Slurm.

1. **No-op (shared FS).** Today. Transport is nothing; the remote path *is* the local path.
2. **Explicit staging (copy in/out).** A `Transport` peer to `Executor` — `stage_in` before submit,
   `stage_out` after. No-op on a shared FS; copy/upload otherwise. Pulsar (`transfer`/`copy`/`none`)
   and Nextflow's FS abstraction are this. Additive; the materializer, ingest, registry, provenance,
   and records do not move.
3. **Fake a POSIX FS over object storage (Fusion-style FUSE).** A thin client mounts S3/etc. as
   POSIX with lazy background download/upload — the worker *thinks* it has a shared FS. Seqera
   benchmarks it "on par with Lustre at the cost of object storage." Note it works by *faking* the
   shared FS, because moving the data is the thing you want to avoid — infrastructure cost, zero
   worker change.
4. **Content-addressed staging (Bazel Remote Execution / CAS).** Inputs and outputs live in a
   content-addressable store keyed by digest; the worker fetches inputs by hash, runs hermetically,
   pushes outputs by hash — "workers don't need mounted shared storage or NFS." Buys dedup and
   *incremental upload* (`FindMissingBlobs` ships only novel blobs — the way to stop re-sending the
   same matrix). **Foothold:** MDV already computes a `content_hash` for analysis identity (ADR-0006);
   CAS uses content digests as both cache key *and* transport address, so that hashing is the seed of
   this path if it is ever taken.

The one option that is a genuinely **different architecture, not a seam**, is a **pull/agent model**
(Globus Compute, HTCondor, a Celery worker living on the cluster): a long-lived agent on the cluster
pulls jobs from the owner. It trades the sealed-courier worker for a connected agent. Even there the
*task function* stays sealed and everything above the executor is untouched, so it is a **sibling
backend**, not a rewrite — but it is the case that does not reduce to a `Transport` swap, and it
earns its own ADR if "truly remote, unmountable" becomes a first-class target.

## Considered and rejected (for now)

- **Make the base `SlurmExecutor` handle both shared-FS and staging.** Bloats the thin CLI class and
  breaks the courier seal by pushing transfer logic toward the worker. Rejected: keep the base
  executor assuming a shared FS; a staging/REST executor is a *separate* backend that owns the
  transfer subsystem.
- **Object storage (S3/MinIO) as "the shared FS."** Not POSIX — you cannot `--chdir` into a bucket
  or `stat` a marker — so the worker would have to pull-from/push-to the bucket itself, breaking the
  "worker only touches its local workspace path" seal. It is network staging wearing a filesystem
  costume; if you want the POSIX illusion, that is option 3 (FUSE), done properly.
- **SSH-remote submit as the remote answer.** Running `sbatch` over SSH from an off-cluster owner
  works but is what the ecosystem is moving *away* from (Parsl is removing its SSH channels in favor
  of co-location / Globus Compute). It also does not solve data movement — you still need a shared FS
  or staging underneath. The injected command-runner already makes local-vs-SSH a `run`-callable
  swap, so nothing is lost by not designing around SSH now.

## Consequences

- The first `SlurmExecutor` stays a thin CLI class that moves no bytes; **"shared POSIX FS at a
  matching path" is its stated precondition**, satisfied by deployment (same node, NFS/Lustre mount,
  managed cloud FS), not by code.
- The `Transport`/staging seam remains **named but unbuilt** (as ADR-0004 already anticipated); when
  a true no-shared-FS deployment arrives it is added as options 2–4 above **without touching the
  executor contract**.
- The pull/agent model is recorded as a **sibling backend** for the unmountable-remote case, with its
  own future ADR — not a reason to change the current design.
- `content_hash` is the foothold for a future CAS transport, should the ladder reach option 4.

Prior art:
[Bazel Remote Execution API (CAS, sealed action)](https://github.com/bazelbuild/remote-apis),
[Seqera — storage architecture for Nextflow pipelines](https://seqera.io/blog/selecting-the-right-storage-architecture-for-your-nextflow-pipelines/),
[Nextflow Fusion file system](https://nextflow.io/docs/latest/fusion.html),
[Pulsar staging (`transfer`/`copy`/`none`)](https://pulsar.readthedocs.io/en/latest/galaxy_conf.html),
[Parsl — removing SSH channels](https://parsl-project.org/2024/10/29/removing-channels.html),
[HPC scratch staging best practice (UIowa)](https://uiowa.atlassian.net/wiki/spaces/hpcdocs/pages/76513434/Scratch+Filesystems),
[Data-locality-aware task scheduling (arXiv)](https://arxiv.org/html/2407.08584v2).
