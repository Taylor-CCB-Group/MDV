# CHATMDV Product Requirements Document (PRD)

## Purpose

ChatMDV enables users to ask natural-language questions about MDV datasets and receive:

- analysis outputs in chat (text and bounded tables), and
- persisted visualizations/charts in MDV when appropriate.

This PRD defines the architecture boundaries, responsibilities, and quality requirements for ChatMDV behavior.

## Scope and Boundaries

ChatMDV behavior is owned by code in `python/mdvtools/llm/` plus selected policy/regression tests.

### In Scope

- Prompting and generation policy
- Datasource-aware field guidance
- Code execution wrapping and diagnostics
- Verification and user-facing validation summaries
- LLM helper utilities used by the orchestration flow

### Out of Scope (Without Explicit Product Decision)

- Core MDV runtime changes in `src/datastore/`
- `mdvproject` data model changes
- Chart manager/runtime infrastructure changes

When a chat-generated chart fails because of bad chart `param` values, first fix prompt/field-discovery policy in ChatMDV components rather than changing core datastore behavior.

## System Components

### 1. Orchestration

**Primary file:** `python/mdvtools/llm/langchain_mdv.py`

**Provider selection:** `python/mdvtools/llm/llm_providers.py` — discovers OpenAI and Ollama models, pairs chat models with embedding backends for RAG, and constructs LangChain LLM/embedding instances.

Responsibilities:

- Coordinate prompt construction and LLM execution
- Run the multi-step pipeline: **pandas agent → RAG codegen → preflight → execute**
- Validate agent and RAG output before execution (no silent stub success)
- Route generated code into execution wrapper
- Attach execution diagnostics back to chat flow
- Keep orchestration logic separate from policy text and verification copy

#### Chat request pipeline (orchestration blocks)

| Block | Role |
|-------|------|
| `b8` / `b10b` | Pandas agent inspects `df1`/`df2`; emits `fields` / `charts` plan |
| `b11` / `b12` | RAG retrieval + codegen LLM |
| `b13` / `b13b` | `prepare_code`, preflight validate (single retry) |
| `b14` | Execute generated script |

Each request logs `model_id`, `provider`, chat model name, and embedding model to `chat_debug.log` and the per-project chat log.

**Agent plan guards** (all providers):

- Resolve `fields` / `charts` from agent output or ReAct intermediate steps (`_resolve_agent_plan`).
- If still missing and the user question is non-empty: **warn and continue** with a question-derived fallback plan (`_fallback_agent_plan`) so RAG retrieval can proceed.
- If the question is empty: fail with a clear error.

**Codegen guards** (all providers):

- Empty RAG `result`: one retry with an explicit fenced-Python instruction, then fail.
- No extracted Python / warning stub from `prepare_code`: fail before `execute_code`.

**Ollama-specific:**

- Functions agent first; ReAct fallback on tool errors or missing field/chart plan.
- Compact RAG prompt (`get_createproject_prompt_RAG(..., compact=True)`) — shorter than full prompt but retains project context, datasource roles, do-not-recreate-datasource policy, and output-format rules.

**OpenAI-specific:**

- Functions agent only (no ReAct fallback).
- Full RAG prompt (`compact=False`).

### 2. Prompt and RAG Policy

**Primary file:** `python/mdvtools/llm/templates.py`

Responsibilities:

- Define generation instructions and retrieval-augmented policy text
- Encode high-level constraints for analysis vs visualization behavior
- Ensure model defaults favor safe, bounded outputs in chat when chart persistence is not valid

#### RAG code corpus

Prompt text lives in `templates.py`; **retrieval** additionally indexes a **local code corpus** of sample MDV scripts so embeddings can surface chart-construction patterns consistent with AnnData-style projects.

- **Location:** `python/mdvtools/test_projects/RAG_examples/`
  - `ANNDATA_examples/` — default indexed corpus (chart recipes for `.h5ad` / AnnData workflows).
  - `PBMC3K_examples/`, `TABULAR_examples/`, `TAURUS_examples/` — extra runnable demos; excluded from default ChatMDV indexing unless callers pass another `directory_path` into the crawler.

- **Indexing:** `python/mdvtools/llm/local_files_utils.py` — `crawl_local_repo()` defaults to `python/mdvtools/test_projects/RAG_examples/ANNDATA_examples/`.

- **When loaded:** Importing `langchain_mdv.py` triggers crawling, chunking, and FAISS index construction on startup.

- **Automation:** `python/mdvtools/test_projects/llm_automated_testing.py` uses the same default `crawl_local_repo()`. `llm_automated_testing_command_line.py` points at the same `ANNDATA_examples` folder but duplicates crawl logic with a different ignore list, so indexed files may not match orchestration exactly.

- **Operational impact:** Adding or changing scripts under `ANNDATA_examples/` alters retrieval behavior after process restart / re-embed. Sphinx API generation excludes `**/RAG_examples/**` (`docs/maindocs/conf.py`).

### 3. Datasource and Field Policy

**Primary file:** `python/mdvtools/llm/datasource_roles.py`

Responsibilities:

- Define datasource-role specific chart/field rules
- Prevent cross-datasource field-id misuse
- Express policy for marker workflows and visualization consistency

Key policy entry points:

- `format_marker_gene_scanpy_fallback_policy()`
- `format_obs_table_chart_param_policy()`
- `format_visualization_consistency_policy()`
- `format_marker_ranking_viz_policy()`
- `CHAT_RANK_GENES_DATASOURCE_NAME`

### 4. Execution Diagnostics

**Primary files:** `python/mdvtools/llm/code_execution.py`, orchestration hooks in `langchain_mdv.py`

Responsibilities:

- Execute generated analysis code in controlled wrapper(s)
- Provide actionable diagnostics for failures
- Differentiate policy failures from runtime failures

### 5. Verification Copy and User Validation

**Primary file:** `python/mdvtools/llm/verification.py`

Responsibilities:

- Produce post-analysis verification summaries for user trust
- Describe what was computed, persisted, or printed
- Guide users to validate that chart outputs match analysis claims

### 6. Supporting Helpers

**Location:** other modules in `python/mdvtools/llm/` (for example `code_manipulation.py`, `column_field_resolve.py`)

Responsibilities:

- Utility transformations and resolver logic
- Optional post-save validation, such as dropping charts whose unresolved `param` tokens do not map to datasource field IDs

### 7. Runnable CLI Module

**Primary file:** `python/mdvtools/llm/chat_cli.py`

Responsibilities:

- Provide one-shot ChatMDV execution via module invocation (`python -m mdvtools.llm.chat_cli`)
- Accept required request inputs (`--project`, `--sentence`)
- Support optional operational controls (`--json`, `--verbose`, `--output-dir`, `--view-name`)
- Return a stable result contract for terminal and automation use
- Persist debug artifacts (`generated_code.py`, `result.json`, `view_snapshot.json`) when `--output-dir` is provided
- Support standalone CLI execution by installing a no-op SocketIO shim when app SocketIO is not initialized

## End-to-End Request Flow

1. User asks a question in chat (optionally selects model from discovered OpenAI/Ollama list).
2. Orchestration resolves chat + embedding models; logs provider/model ids.
3. Pandas agent inspects datasources and produces a field/chart plan (with Ollama ReAct fallback and question-derived fallback when needed).
4. RAG prompt is built (compact for Ollama, full for OpenAI); retrieval augments with `ANNDATA_examples` corpus.
5. Codegen LLM returns fenced Python; guards reject empty output or stubs.
6. `prepare_code` wraps script; preflight validates chart API usage (single retry on failure).
7. Execution wrapper runs code and captures diagnostics.
8. Chat response includes text/tables and optional persisted chart instructions.
9. Verification layer summarizes what to check for correctness and consistency.
10. Optional field/token resolution logic validates persisted chart parameters.

### Runnable Module Flow (`python -m`)

1. User runs `python -m mdvtools.llm.chat_cli --project <path> --sentence "<prompt>"`.
2. CLI validates project path and datasource availability.
3. CLI invokes `ProjectChat.ask_question(...)` and receives `AskQuestionResult`.
4. If `--view-name` is provided, CLI renames the generated view to the requested name.
5. CLI returns success/failure with:
   - `success`
   - `project_path`
   - `views_file`
   - `view_name`
   - `message`
   - `debug_output_dir`
6. If `--output-dir` is provided, CLI writes debug artifacts for inspection.

## Critical Policy Rules

### Marker Genes with Missing DE Columns

If `.h5ad` exists at `data_path`, model should compute markers via Scanpy and print bounded results in chat instead of assuming marker columns already exist on `genes`.

### Table Charts vs Chat-Only Tables on `cells`

Chart `param` tokens must be valid field IDs on the chart's target datasource.  
For top-marker Scanpy outputs, default to bounded `print(...)` unless required columns already exist on the target datasource.

### Visualization/Analysis Consistency

Persisted charts must reflect the same computational pipeline as the printed analysis output. Avoid mixing incompatible pipelines for the same quantitative claim.

### `rank_genes_groups` and Observation Datasource

Do not substitute wrapper-based DotPlot/Heatmap on `cells` for full Scanpy marker statistics tables. Prefer bounded print output or persist long-format results with:

- `add_datasource('chat_rank_genes_result', ..., replace_data=True, add_to_view=...)`

Avoid `set_view` patterns that overwrite `initialCharts` and drop scratch analysis tables.

## Test Strategy

Policy and regression tests should be concentrated in:

- `python/mdvtools/tests/test_chat_first_text_table_policy.py`
- `python/mdvtools/tests/test_chat_ollama_pipeline.py` — agent plan resolution, ReAct fallback, RAG empty guard, stub rejection, compact prompt
- `python/mdvtools/tests/test_llm_providers.py` — provider discovery, Ollama URL defaults, embedding pairing
- `python/mdvtools/tests/test_langchain_mdv_helpers.py` — text/table intent helpers, init failure when no providers
- `python/mdvtools/tests/test_datasource_roles.py`
- `python/mdvtools/tests/test_column_field_resolve.py`
- `python/mdvtools/tests/test_code_execution.py`
- `python/mdvtools/test_projects/chatmdv_csv_eval/run_prompt_csv.py` — pbmc3k CSV benchmark harness (10 prompts, subprocess row isolation)
- `python/mdvtools/test_projects/llm_automated_testing.py` — end-to-end LLM harness using the same default RAG corpus as orchestration (`RAG_examples/ANNDATA_examples/`).

Add new ChatMDV policy tests in these suites rather than unrelated runtime test areas.

## Change Management and Drift Control

Before merge, validate boundary compliance:

```bash
python python/mdvtools/llm/check_chatmdv_boundary.py
python python/mdvtools/llm/check_chatmdv_boundary.py --git-range origin/main...
```

A non-zero exit indicates files changed outside approved ChatMDV boundaries.

## Debugging Workflow

1. Localize failure: prompt policy vs execution diagnostics vs verification copy.
2. Change the smallest responsible file in the ChatMDV boundary.
3. Add or extend a focused policy/regression test.
4. Run targeted tests; run boundary check if scope drift is possible.

## Non-Functional Requirements

- **Predictability:** deterministic policy constraints over ad hoc behavior.
- **Safety:** prefer bounded printed outputs when persistence preconditions are unmet.
- **Traceability:** diagnostics and verification copy should make failures understandable.
- **Separation of concerns:** policy, execution, and verification remain independently evolvable.

## Definition of Done (for ChatMDV Changes)

- Change remains within ChatMDV boundary (or has explicit product approval to cross it).
- Relevant policy/regression tests are added or updated.
- Generated behavior follows datasource-field constraints and marker workflow policy.
- Verification text accurately reflects what was executed and persisted.

## Local Ollama and multi-provider LLM

ChatMDV supports **OpenAI** (when `OPENAI_API_KEY` is set) and **Ollama** (OpenAI-compatible `/v1` API). The Chat UI exposes a model dropdown; CLI and CSV eval accept `--model <id>` or `CHATMDV_DEFAULT_MODEL`.

### Recommended models

| Role | OpenAI | Ollama |
|------|--------|--------|
| Chat / codegen | `gpt-4.1`, `gpt-4o` | `qwen2.5-coder:7b`, `qwen2.5-coder:14b`, `qwen3-coder:latest` |
| Embeddings (RAG) | `text-embedding-3-large` | `nomic-embed-text` (`ollama pull nomic-embed-text`) |

Avoid general chat models (e.g. `gemma4:26b`) for the codegen pipeline; they often return empty agent/RAG output or hit ReAct iteration limits.

### Environment

| Variable | Purpose |
|----------|---------|
| `OPENAI_API_KEY` | Enables OpenAI chat + embeddings |
| `OLLAMA_BASE_URL` | Ollama host (Docker devcontainer default: `http://host.docker.internal:11434`) |
| `CHATMDV_DEFAULT_MODEL` | Optional default chat model id (e.g. `ollama:chat:qwen3-coder:latest`) |
| `CHATMDV_CSV_EVAL_OUTPUT_DIR` | Default directory for CSV benchmark artifacts |

Verify Ollama from the container:

```bash
curl -s http://host.docker.internal:11434/api/tags
```

### Ollama pipeline behavior

When the selected model provider is Ollama:

1. **Functions agent first** — inspects `df1`/`df2` via tool calling (~3–5s typical).
2. **ReAct fallback** — only on functions-agent exception or missing `fields`/`charts` plan.
3. **Agent plan resolution** — extract plan from ReAct intermediate steps when final output hits iteration limits.
4. **Question fallback** — if both agents fail but the user question is non-empty, log a warning and continue with `_fallback_agent_plan(question)` (chart hints inferred from question keywords).
5. **Compact RAG prompt** — shorter prompt for local models; includes do-not-`add_datasource` policy for existing project datasources.
6. **Fail fast** — empty RAG result (after one retry) or missing extracted Python blocks execution (no silent `# WARNING:::: No code captured` success).

OpenAI uses the same agent/codegen guards but **no ReAct fallback** and the **full** RAG prompt.

### Alignment expectations

| In scope (current pipeline) | Out of scope (later work) |
|-----------------------------|---------------------------|
| Reliable codegen or clear errors | Curated RAG corpus expansion |
| Field/chart plan from agent REPL + fallbacks | MDV-native tools / MCP |
| Compact prompts + guards for local models | Fine-tuning on project Q→code traces |
| Fail-fast on empty LLM output / stubs | vLLM or non-Ollama local serving |

### CSV eval benchmarks

Track pipeline quality with the pbmc3k prompt CSV harness. Project: `/app/mdv/pbmc3k_chat` (create via `setup_pbmc3k_chat.py`).

```bash
cd python
export CHATMDV_CSV_EVAL_OUTPUT_DIR=/app/logs/chatmdv
mkdir -p "$CHATMDV_CSV_EVAL_OUTPUT_DIR"

# Smoke (one row)
.venv/bin/python mdvtools/test_projects/chatmdv_csv_eval/run_prompt_csv.py \
  --csv mdvtools/test_projects/chatmdv_csv_eval/prompts_pbmc3k.csv \
  --model qwen3-coder:latest \
  --output-dir "$CHATMDV_CSV_EVAL_OUTPUT_DIR" \
  --limit 1

# Full run (Ollama)
.venv/bin/python mdvtools/test_projects/chatmdv_csv_eval/run_prompt_csv.py \
  --csv mdvtools/test_projects/chatmdv_csv_eval/prompts_pbmc3k.csv \
  --model qwen3-coder:latest \
  --output-dir "$CHATMDV_CSV_EVAL_OUTPUT_DIR"

# Full run (OpenAI baseline)
.venv/bin/python mdvtools/test_projects/chatmdv_csv_eval/run_prompt_csv.py \
  --csv mdvtools/test_projects/chatmdv_csv_eval/prompts_pbmc3k.csv \
  --model openai:chat:gpt-4.1 \
  --output-dir "$CHATMDV_CSV_EVAL_OUTPUT_DIR"
```

**Artifacts** (timestamped under `CHATMDV_CSV_EVAL_OUTPUT_DIR`):

| File | Use |
|------|-----|
| `benchmark_summary_*.json` | `exec_success_rate`, `view_has_charts_rate`, `by_complexity` |
| `timing_*.jsonl` | Per-row block timings (`block_b10b_s`, `block_b12_s`, `block_b14_s`) |
| `results_*.csv` | Row-level compact metrics |
| `failures_*.jsonl` | Failed rows with `failure_reason` |

**Reference results** (pbmc3k, 10 prompts, 2026-06-13):

| Config | `exec_success_rate` | `view_has_charts_rate` | Notes |
|--------|---------------------|------------------------|-------|
| Ollama `qwen3-coder:latest` — pre-guard baseline | 80% | 80% | Functions agent, full prompt, no stub guards |
| Ollama — ReAct default (regression) | 30% | 30% | ReAct-first caused iteration-limit failures |
| Ollama — functions + compact + strict agent guard | 50% | 50% | Fail-fast on missing plan hurt success rate |
| Ollama — **current** (functions + compact + question fallback + guards) | **90%** | **80%** | Run `20260613T102009Z`; 1 preflight failure (row 9, PCA) |
| OpenAI `gpt-4.1` | 90% | 90% | Run `20260613T084756Z`; 1 runtime failure (row 10, GNLY) |

**Design notes from benchmarking:**

- ReAct as the **default** Ollama agent regressed success (80% → 30%); functions-first with ReAct **fallback** is required.
- Strict fail-fast on missing agent plan regressed Ollama (80% → 50%); question-derived fallback restores parity while keeping stub/RAG guards.
- Compact prompt alone did not beat full prompt; compact **with** do-not-recreate-datasource policy is the chosen Ollama default.
- Remaining Ollama failures are mostly **codegen/preflight** quality (invalid chart kwargs), not pipeline guard false positives.

## Runnable Module Specification

### Invocation

```bash
python -m mdvtools.llm.chat_cli --project <project_path> --sentence "<prompt>"
```

### Arguments

- Required:
  - `--project`: MDV project directory
  - `--sentence`: one-shot natural-language request
- Optional:
  - `--output-dir`: write debug artifacts to directory
  - `--json`: print machine-readable output
  - `--verbose`: print additional human-readable details
  - `--view-name`: force final persisted view name

### `--view-name` Semantics

- If provided and different from the auto-generated view name:
  - copy generated view payload to the requested name
  - delete original auto-generated view entry
  - return requested name in output `view_name`
- Empty/whitespace-only value is invalid.

### Output Contract

```json
{
  "success": true,
  "project_path": "/abs/path/to/project",
  "views_file": "/abs/path/to/project/views.json",
  "view_name": "custom_or_generated_view_name",
  "message": "Success",
  "debug_output_dir": "/tmp/chatmdv-debug-run1"
}
```
