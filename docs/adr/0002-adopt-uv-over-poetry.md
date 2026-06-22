# Adopt uv (build backend + lockfile) over Poetry

**Status:** accepted — landed. uv is on `main` via PR #472; this branch consolidates on top of it.

## Context & decision

Consolidating to a single `pyproject.toml` (ADR-0001) requires expressing
`[project.optional-dependencies]`, which needs **static PEP 621 `[project]` metadata** —
something Poetry's `dynamic = [...]` stubs + legacy `[tool.poetry.*]` groups could not do.

A uv migration was attempted (**PR #468**), reverted (**PR #471**, no recorded rationale),
then **re-landed cleanly on `main` as PR #472** while this work was in progress. #472
converted the manifest to static PEP 621, switched the build backend to `uv_build`, added
`uv.lock`, pinned uv (`[tool.uv] required-version`), and moved the Dockerfile/CI to uv.

We **adopt uv** (build backend `uv_build` + lockfile/installer `uv`; Poetry dropped) and
build the packaging consolidation *on top of* #472 — rather than re-doing the migration
ourselves.

## Why

- The consolidation and the uv migration are the **same refactor** — both hinge on the
  PEP 621 conversion. #472 did that conversion; we only add the `app` extra + guards.
- It matches the team's stated direction (off Poetry).
- Path 1 (ADR-0001) means we never needed Hatchling's surgical `force-include`; `uv_build`'s
  basic include/exclude (notebooks, scratch JSON, test data, `*.map`) is sufficient.

## Considered and rejected

- **Keep Poetry** (Poetry 2.x supports PEP 621 extras; ~2-line Dockerfile change). The
  *original proposal's* recommendation, chosen to avoid reopening the uv question. Rejected:
  leaves us on Poetry against the team's preference, and #472 settled the question anyway.
- **Consolidate on Poetry now, switch to uv later.** Touches `pyproject.toml` twice. Rejected.

## Consequences

- The earlier "**owed:** why was #471 reverted?" risk is now **moot** — #472 re-landed uv on
  `main` and is passing CI, so we are no longer re-landing blind. (If the original revert
  reason still matters historically, Peter Todd has it; it no longer gates this work.)
- #472 modelled the optional stack as PEP 735 `[dependency-groups]` (`dev`/`docs`/`backend`/
  `llm`/`auth`) and **kept the chat stack + `ruff` in core `dependencies`** — i.e. a bare
  `pip install mdvtools` from `main` was *not* slim. This branch finishes the job: the
  server/chat/auth runtime deps move into a single `[project.optional-dependencies].app`
  extra; `backend`/`llm`/`auth` groups are deleted; `dev`/`docs` groups remain (dev-only).
- Dockerfile install lines change from `uv sync … --group dev --group backend --group auth`
  to `uv sync … --group dev --extra app` (both the cache layer and the final layer).
- Editing the dependency layout invalidates `uv.lock`; it must be regenerated with `uv lock`
  (in the devcontainer, with uv at the pinned version). The regenerated lock showed **no
  version or package changes** — purely re-classification.
- **`ruff` is a runtime dependency, not dev-only** (see ADR-0001 / proposal §4): the chat
  feature shells out to `ruff check` to lint generated code, so it lives in the `app` extra
  **and** the `dev` group (the latter for `make lint`/`make format`).
- The `[tool.uv] required-version` is an **exact** pin (`==`), which makes every small uv
  drift fail `uv lock`/`uv sync`. Flagged as a follow-up (consider `>=x,<y`); not changed
  here since it came in with #472.
