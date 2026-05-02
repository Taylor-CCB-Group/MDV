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

Responsibilities:

- Coordinate prompt construction and LLM execution
- Route generated code into execution wrapper
- Attach execution diagnostics back to chat flow
- Keep orchestration logic separate from policy text and verification copy

### 2. Prompt and RAG Policy

**Primary file:** `python/mdvtools/llm/templates.py`

Responsibilities:

- Define generation instructions and retrieval-augmented policy text
- Encode high-level constraints for analysis vs visualization behavior
- Ensure model defaults favor safe, bounded outputs in chat when chart persistence is not valid

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

1. User asks a question in chat.
2. Orchestration builds context and applies template/policy guidance.
3. LLM generates analysis code and output directives (print/chart intent).
4. Execution wrapper runs code and captures diagnostics.
5. Chat response includes text/tables and optional persisted chart instructions.
6. Verification layer summarizes what to check for correctness and consistency.
7. Optional field/token resolution logic validates persisted chart parameters.

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
- `python/mdvtools/tests/test_datasource_roles.py`
- `python/mdvtools/tests/test_column_field_resolve.py`
- `python/mdvtools/tests/test_code_execution.py`

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
