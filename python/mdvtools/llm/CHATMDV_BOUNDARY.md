# ChatMDV boundary

Chat-generated behavior (prompting, code policy, execution wrappers, verification text) lives under **`python/mdvtools/llm/`** and selected tests. **Do not change core MDV runtime** (`src/datastore/`, `mdvproject` data model, chart manager) to fix Chat issues unless there is a separate, explicit decision.

## Allowed areas (default)

| Area | Paths |
|------|--------|
| Orchestration | [`langchain_mdv.py`](langchain_mdv.py) |
| Prompts / RAG | [`templates.py`](templates.py) |
| Datasource hints for LLM | [`datasource_roles.py`](datasource_roles.py) |
| Subprocess execution wrapper | [`code_execution.py`](code_execution.py) |
| Chat verification summaries | [`verification.py`](verification.py) |
| Supporting LLM helpers | Other modules in this directory (e.g. `code_manipulation.py`, `column_field_resolve.py`) |

## Tests (regression / policy)

- [`python/mdvtools/tests/test_chat_first_text_table_policy.py`](../tests/test_chat_first_text_table_policy.py)
- [`python/mdvtools/tests/test_datasource_roles.py`](../tests/test_datasource_roles.py)
- [`python/mdvtools/tests/test_code_execution.py`](../tests/test_code_execution.py)

Add new ChatMDV policy tests here rather than spreading across unrelated test files.

## Separation of concerns

1. **Policy** (what the model should generate): `templates.py`, `datasource_roles.py`.
2. **Execution diagnostics** (how failures are reported): `code_execution.py`, `langchain_mdv.py` (orchestration only).
3. **Verification copy** (what the user sees to validate results): `verification.py`.

If the UI errors on a bad chart `param` (e.g. wrong field id), **fix prompts / field discovery first**; do not alter the datastore unless the product owner approves infrastructure changes.

**Marker genes / missing DE columns:** Use `format_marker_gene_scanpy_fallback_policy()` in [`datasource_roles.py`](datasource_roles.py) and the RAG text in [`templates.py`](templates.py). When a `.h5ad` exists at `data_path`, the model should compute markers in Scanpy and print to chat instead of asserting cell-level columns on the `genes` datasource.

**Visualization vs analysis consistency:** Use `format_visualization_consistency_policy()` in [`datasource_roles.py`](datasource_roles.py). Saved MDV charts must reflect the same pipeline as printed tables—do not mix Scanpy/AnnData outputs with unrelated wrapper-based expression heatmaps for the same quantitative claim.

## Pre-merge drift check

From repo root:

```bash
python python/mdvtools/llm/check_chatmdv_boundary.py
```

With a base ref (e.g. before opening a PR):

```bash
python python/mdvtools/llm/check_chatmdv_boundary.py --git-range origin/main...
```

Exits with code `1` if any changed file falls outside the boundary lists.

## Debugging workflow (short)

1. Localize: prompt vs execution vs verification.
2. Change the smallest file in the table above.
3. Add or extend a test under `python/mdvtools/tests/test_chat_*` or `test_datasource_roles.py`.
4. Run targeted pytest, then run the boundary script if you touched non-`llm/` files accidentally.
