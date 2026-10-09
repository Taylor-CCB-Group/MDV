# Resource requirements and node placement: matching a job's needs to a node's capabilities

**Status:** proposed — design captured, **build deferred (parked)**. No consumer needs it yet;
`SlurmExecutor` currently emits basic directives only. This ADR records the seam so the design is
settled before the first tool (a large UMAP, or the GPU backend) forces it.

## Context & decision

Slice A of `SlurmExecutor` deliberately renders **basic directives only** — `--job-name`,
`--chdir`, `--output`/`--error` — and nothing about *how much machine* the job needs. Real jobs
differ: a UMAP over a large matrix wants more memory; the parked GPU backend
(rapids-singlecell, ADR-0008/[[project_mdvtools_packaging]]) needs an actual GPU node. Expressing
that is a **new framework primitive**: matching a job's resource **requirement** (mem / cpus /
walltime / gpus) to a node **capability**. It is genuinely new surface — not a variation on
anything the courier or dispatch ADRs already cover — so it gets its own record.

The decision this ADR fixes: **a job declares abstract resource requirements; the executor
translates them into whatever its scheduler speaks.** `--gres=gpu:1` on Slurm,
`resources.limits: nvidia.com/gpu: 1` on Kubernetes, `num_gpus=1` on Ray — same requirement,
per-backend dialect. The framework carries the requirement in neutral terms and never hard-codes a
scheduler flag above the executor.

## Where a requirement comes from

Two sources, layered:

- **Per-tool default on `ToolSpec`** — a `resources` field (e.g. `{"mem": "8G", "cpus": 4,
  "time": "01:00:00", "gpus": 0}`). This is the tool author saying "UMAP-over-a-matrix generally
  needs this much." It travels with the spec, like `input_shape` and `entrypoint` already do, and
  is the registry's answer to placement just as it is already the dispatch table.
- **Per-submission override (optional)** — a caller may raise the ask for a known-large input
  (this dataset is 10× the usual). Overrides the tool default for that one job.

Provisional home: `ToolSpec.resources` as the base, per-submission override on top. Not locked —
there is no consumer yet, so the exact shape is finalized when the first tool actually needs
non-default resources. The GPU story is the likely first mover.

## How the requirement reaches the executor

Today the seam is `Executor.submit(entrypoint, workspace)`. Resources have to travel
spec → manager → `submit`, so the **Protocol widens** to:

```
submit(self, entrypoint: str, workspace: Path, resources: dict | None = None) -> Handle
```

- `LocalSubprocessExecutor` **ignores** `resources` — a bare subprocess has no placement; whatever
  the box has is what it gets. (This is also why the parameter is optional-with-default: Local's
  behavior is unchanged and existing callers don't break.)
- `SlurmExecutor` turns `resources` into sbatch directives in `_render_script`: `mem → --mem`,
  `cpus → --cpus-per-task`, `time → --time`, `gpus → --gres=gpu:N`. This is exactly where the
  directives slice A left out finally land.
- A future `K8sExecutor`/`RayExecutor` maps the same dict onto its own request/limit vocabulary.

The neutral `resources` dict is the contract; the per-scheduler translation is localized in each
executor — the same shape as `poll`'s running/done/lost translation (ADR-0008): one abstract
vocabulary, N backend dialects, no cross-backend branching in the manager.

## This folds in the GPU backend

`--gres=gpu:N` is not a separate feature from GPU support — it *is* how a GPU-backed worker
requests a GPU node. So the parked GPU work splits cleanly along this seam:

- **Placement half** (this ADR): `gpus: 1` in `resources` → `--gres=gpu:1`. Framework concern.
- **Compute half** (parked with the worker): a `backend` flag selecting cuML/rapids-singlecell vs
  umap-learn, the lazy optional `gpu` extra. Tool concern.

⚠️ The `backend` flag **must feed the content-hash** (ADR-0006): cuML and umap-learn produce
*different coordinates* for the same input, so a run on GPU is a different analysis identity from a
run on CPU. Placement (which node) and identity (which algorithm) are separate — only the latter
enters the hash.

## What is built vs deferred

- **Built:** basic directives only (slice A). No `resources` parameter on `submit` yet.
- **Deferred (this ADR):** the `resources` field, the `submit` widening, the `_render_script`
  directive mapping, and the GPU placement path. The testable slice when built mirrors slice A's
  directive tests — feed a `resources` dict, assert `--mem`/`--cpus-per-task`/`--time`/`--gres`
  appear in the rendered script.

Recorded now, not built now, because nothing consumes it yet and the first real consumer (large
UMAP or the GPU backend) will pin the exact `resources` shape better than guessing ahead of it.

## Consequences

- The `submit` Protocol gains an optional `resources` param; Local ignores it, so nothing regresses
  before the Slurm mapping is built.
- Resource declaration lives on the tool (with an optional per-submission override), keeping the
  registry the single place that describes a tool's dispatch *and* placement.
- The GPU backend stops being one monolithic parked item and becomes two: a placement path that is
  just another entry in the `resources` dict, and a compute-backend flag that is a worker/hash
  concern.

Prior art: [Slurm `--gres` / consumable resources](https://slurm.schedmd.com/gres.html),
[Kubernetes resource requests & limits](https://kubernetes.io/docs/concepts/configuration/manage-resources-containers/),
[Nextflow process resources (`cpus`/`memory`/`accelerator`)](https://www.nextflow.io/docs/latest/process.html#resources),
[Ray resource requirements (`num_gpus`)](https://docs.ray.io/en/latest/ray-core/scheduling/resources.html).
