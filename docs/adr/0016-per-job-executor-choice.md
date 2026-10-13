# Per-job executor choice: the user picks where a job runs, the operator picks how Slurm is reached

**Status:** accepted (jobs framework).

## Context

ADR-0008 made executors peers behind one seam (`submit`, `poll`, `locate_result`), and the
`SlurmExecutor` was validated on a real cluster. Until now each project's `JobManager` held one
executor, chosen in code, and the server always built managers with `LocalSubprocessExecutor`. A
deployment such as BMRC, Oxford's HPC, can run a small job on the server and send a large one to
Slurm, so the person submitting the job needs to choose. BMRC also offers the Slurm REST API
(`slurmrestd`) as well as the CLI, so the framework needs a second way to reach Slurm.

## Decision

The user chooses where a job runs, from the executors the deployment offers, when they submit it.
The Run analysis dialog shows a dropdown only when the server offers more than one executor, and
then the job list also shows a "Run on" column with each job's executor, where a record with no
executor shows the default. The
operator decides once per deployment how the server reaches Slurm (the CLI with `sbatch` and
`squeue`, or `slurmrestd`), because that depends on what the cluster allows and the user cannot act
on it.

A `JobManager` holds a named set of executors. Each entry carries its own executor, workspace root
and concurrency limit, because local scratch and the shared Slurm scratch (ADR-0010) are different
directories and only the local executor needs a cap. The first entry is the default.

The record stores the executor's name from the moment it is queued (`JobRecord.executor`). A
queued record has no handle yet, so the handle's `kind` cannot carry the choice. The manager
resolves the default at submit and stores the name, so a later change to the default does not move
jobs already queued. Records written before this field existed load with `executor = None` and run
on the default.

Every step that touches a record uses that record's executor: the workspace path, dispatch, the
poll in `tick`, and the poll in the restart reconcile (ADR-0005). Concurrency is counted per
executor, so a full local pool does not hold back a Slurm job.

## Configuration

The server reads the deployment's executors from the environment when it builds a manager:

- `local` is always offered, with the default workspace root and a limit of 2 concurrent jobs
- `MDV_JOBS_SLURM` set to `cli` adds `slurm`, using `SlurmExecutor` with no concurrency limit,
  since Slurm queues the work itself
- `MDV_JOBS_SLURM_WORKSPACE_ROOT` sets the shared scratch root for `slurm`, and is required when
  `MDV_JOBS_SLURM` is set; each project gets its own directory under it

`MDV_JOBS_SLURM=rest` will select the `slurmrestd` executor once it exists. The user-facing name
stays `slurm` in both cases.

## Consequences

- A submit naming an executor the server does not offer is rejected with a 400 and writes no record.
- A record whose executor is no longer configured, for example after the operator removes Slurm,
  is failed with an error naming the executor. The manager does not guess another executor, because
  a Slurm handle means nothing to the local executor.
- `GET /jobs/executors` lists the names the server offers, and `POST /jobs` accepts an optional
  `executor`. The client view of a record includes `executor`, and the handle stays internal
  (ADR-0012).
