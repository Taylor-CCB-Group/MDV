# One backend-owned, data-only tool spec (params + outputs): GuiSpecType-native, hydrated not serialized, with job-keyed idempotent output

**Status:** accepted — POC design decision (DGE jobs framework). ChatMDV-ready, not
ChatMDV-integrated.

## Context & decision

"Registry-driven tools" is a POC success criterion, so each tool declares itself **once, on
the backend**, as a **serializable, data-only** spec — not a hardcoded form. The spec has
three consumers from one source of truth: the **frontend** hydrates it into live controls,
the **backend** validates submitted params against it, and a **future ChatMDV** reads it as
an LLM function-calling schema.

Full tool contract: `id` / `name` / `description` + **input contract** (needs `cells` + an
expression subgroup) + **params spec** + **output spec** + **worker entrypoint** (the registry
doubles as the dispatch table `tool_id → entrypoint`).

## Params spec — GuiSpecType-native, hydrated not serialized

The spec reuses MDV's existing control vocabulary (`GuiSpecType`: `column`,
`single_category_selection`, `dropdown`, `slider`, …) with `columnType` constraints. We chose
this over JSON Schema because the **frontend consumer ships now** and gets MDV column pickers
+ live datasource binding for free; the LLM consumer is deferred, so its **`GuiSpecType →
JSON Schema` projection** (column→field-id, category→enum, slider→number) is a documented
future seam. (Outside bio, the same "one spec → form + validation + LLM" *is* JSON Schema —
MCP tool definitions are JSON Schema. We project to it later rather than adopt it now.)

This is a **new data-only spec the frontend *hydrates into* a live `GuiSpec`** — **not** a
serialization of `GuiSpec`, which is impossible: `GuiSpec` controls carry `func` callbacks and
`sourceColumn: () => FieldSpec` (functions, `charts.d.ts`). The spec carries no callbacks; the
frontend attaches them on arrival.

**Inter-param bindings: exactly two kinds, no general DSL.**
1. *options-from-another-param's-column* — `target` / `reference` choices are the categories of
   the column named by `groupby`.
2. *visibility-conditional-on-mode* — the comparison-mode selector (column-vs-value modes show
   `groupby`/`target`/`reference`; "use current selection" shows none).

The backend re-validates both at submit (are `target`/`reference` real categories of `groupby`
**in the filtered subset**; are both groups non-empty). The format may grow a third binding
when a second tool actually needs one — not before. (This is Galaxy tool-XML's
`<conditional>` / `data_column` territory, where all the complexity lives; we cap it
deliberately.)

## Safety

The registry is an **allow-list of vetted tools with typed, validated params**. A future LLM
selects a `tool_id` and fills typed blanks; it **never executes free-form code**. Backend
validation is the gate. (Standard function-calling safety — recorded for completeness, not the
novel part.)

## Output spec — a list of typed outputs, identified by the job

The output spec describes a **list of typed outputs**, each declaring its **shape** (gene
table / cell columns) and **target** (new datasource / existing `cells`) — not a single fixed
output. The POC builds the primary output only: a **new standalone datasource per run**, one
row per gene. (Per-cell `obs` write-back is expressible but deferred.)

**Output identity is the `job_id`, fixed at submit — never a timestamp generated at ingest.**
Ingest is **replace-or-skip keyed on `job_id`** (never append). This makes ingest
**idempotent** — the property ADR-0005's recovery depends on: a completion delivered twice
(re-queue *or* at-least-once Socket.IO) yields **one** datasource. A `job_id` also does the
original collision-avoidance job **better** than a timestamp (two submissions in the same
second collide; clock skew/timezones vanish). A human-readable label
(`dge_<groupby>_<target>_vs_<reference>`) is fine for *display*; any timestamp in it must be
the **submit time captured once in the manifest** and reused verbatim — "fixed at submit," not
`now()` at ingest.

## Consequences

- One source of truth eliminates the frontend/backend params drift a separate form + validator
  would invite, and makes the later LLM projection mechanical.
- **Deferred — content-hash result cache.** A key of `hash(tool + params + cells + layer)`
  (Bazel's action-cache idea) enables *"you've already run this exact analysis — open the
  existing result, or run anyway?"* It is **self-invalidating** (editing cells/params changes
  the hash → clean miss) and — because we chose submission-identity — must be an **offer, not a
  silent swap** (preserve "I clicked Submit, I got a result"). Recorded in provenance now; not
  wired up in the POC.

Prior art: [MCP tools (JSON Schema)](https://modelcontextprotocol.io/specification/2025-06-18/server/tools),
Galaxy tool XML (`<conditional>` / `data_column`), [Stripe idempotency keys](https://docs.stripe.com/webhooks).
