# mdvtools Packaging & Release Proposal

**Status:** Decisions locked (grilling session 2026-05-29). See `docs/adr/0001`, `docs/adr/0002`, and `CONTEXT.md`.
**Author:** Marouf Shaikh
**Date:** 2026-05-29
**Goal (set by management):** Stop maintaining two diverging `pyproject.toml` files. Manage **one** well-structured `pyproject.toml` that serves both the slim `pip install mdvtools` (the majority of users) and the full Docker image. Maintainability first; smaller install is a secondary benefit.

**Recommendation (decided):** Collapse to a **single `python/pyproject.toml`** on **uv**. Required deps = slim core; the entire optional cluster becomes **one `app` extra** + guarded imports. `pip install mdvtools` → slim; the Docker image runs `uv sync --all-extras` → full. **Delete `mdvtools_lite_build/`.** No second package, no second version, no `force-include`, ship all modules.

---

## 1. The problem

We maintain **two** `pyproject.toml` files that drift apart. They are not two products — they are two *consumption modes of the same package*:

- **Slim** → `pip install mdvtools`, built today from `mdvtools_lite_build/pyproject.toml` (Hatchling, cherry-picked modules via `force-include`). **The vast majority of users.**
- **Full** → the Docker image, which uses `python/pyproject.toml` (Poetry, all dependency groups). **A small minority.** *Never published to PyPI* — it is only a dependency spec for the image.

Because both files hand-declare `name = "mdvtools"`, `version`, and overlapping dependency lists, every change must be mirrored across both. That redundancy is the divergence:

| Drift point | Today |
|---|---|
| Build backends | Poetry (`python/`) vs Hatchling (`mdvtools_lite_build/`) |
| Dependency lists | two, hand-synced (already disagree on `polars`, `spatialdata`, `pandas`) |
| Version strings | two, both manually `1.2.3` |
| Module selection | a manual `force-include` list reaching into `../python/mdvtools/...` that silently goes stale |

> **Repo history note:** a uv migration (PR #468) was merged then **reverted (PR #471)** with no recorded rationale. The dev side is currently back on Poetry with `dynamic` stubs. We are re-landing uv (§5) — recovering *why #471 reverted* is the one owed item (§7).

---

## 2. Martin's concern, answered with numbers

Martin built the lite wheel and observed:

> *"I use hatch to cherry-pick modules… even doing this, with all the dependencies, the install is still **over 1.4 GB**."*

That observation is the whole argument, and Martin draws the right conclusion himself: *"the only other way is optional dependencies + a graceful way of handling them."* The build backend's file-selection power operates on the axis that doesn't matter (`.py` source), while the 1.4 GB lives in **dependencies**, which `[project.optional-dependencies]` controls — backend-agnostically.

| What | Size | Notes |
|---|---|---|
| `dbutils` + `auth` + `llm` **.py source** (what the lite build actually excludes) | **~576 KB** | the modules in question |
| `mdvtools/static` (built JS frontend) | **59 MB** | already shipped in the slim wheel today, in both builds |
| Full install footprint | **~1.2–1.4 GB** | matches Martin's number |
| └ dominated by deps: faiss, langchain*, scanpy, spatialdata, … | **≈ all of it** | gated by the `app` extra |

**Correction to an earlier framing:** the lite build does **not** exclude `spatial`. Its `force-include` ships `mdvtools/spatial` and it hard-requires `spatialdata`/`spatialdata-io`/`ome-zarr`. So today's PyPI wheel already includes spatial. What the lite build excludes is `dbutils`/`auth`/`llm`. Excluding those saves ~576 KB of ~1.2 GB ≈ **0.05%** — rounding error. (See ADR-0001.)

Every comparable project ships all its modules and guards imports — pandas ships `io.excel` without openpyxl; transformers ships every model without torch; scanpy ships everything. None maintain a second toml to trim KB.

---

## 3. Target design — one `pyproject.toml`, slim by default, on uv

A wheel's **required** deps are the slim core; the single **`app` extra** is pulled only on request. One file serves both modes:

```
pip install mdvtools          → required deps only        = slim  (majority of users)
uv sync --all-extras          → required + app extra      = full  (what the Dockerfile runs)
pip install "mdvtools[app]"   → the pip-equivalent full install
```

```toml
[project]
name = "mdvtools"
version = "1.2.3"                # single source of truth
requires-python = ">=3.11,<3.13"
dependencies = [                 # SLIM CORE — includes spatial + the serve() path
  "anndata==0.12.2", "pandas==2.3.2", "numpy<2.3.0", "numcodecs>=0.16.1",
  "h5py>=3.10.0", "mudata>=0.3.0", "scanpy==1.11.4", "scipy>=1.11.3",
  "Flask==3.0.3", "Flask-SocketIO>=5.3.6", "werkzeug>=3.0.2",
  "polars>=1.35", "fasteners>=0.18", "gevent>=25.5.1",
  # spatial stays in core (preserves current PyPI behaviour):
  "spatialdata>=0.7.2", "spatialdata-io>=0.6.0", "ome-zarr>=0.13.0",
  "generate-tiff-offsets>=0.1.7",   # imported by spatial/jp2k.py — MUST be core
]

[project.optional-dependencies]
# ONE extra: auth+dbutils+llm are a runtime import cycle, so they install together or not at all.
app = [
  "flask-sqlalchemy>=3.1.1", "sqlalchemy", "psycopg2>=2.9.9", "gunicorn==23.0.0", "psycogreen>=1.0.2",
  "langchain>=0.3.27", "langchain-openai>=0.3.33", "langchain-experimental>=0.3.4", "langchain-community",
  "faiss-cpu==1.12.0", "matplotlib>=3.8", "regex", "python-dotenv", "nbformat", "tabulate",
  "gspread>=6.1.2", "oauth2client>=4.1.3",
  "authlib==1.6.4", "python-jose>=3.3.0", "auth0-python>=4.9.0", "flask-jwt-extended>=4.7.1",
  "redis>=5.2.1", "requests>=2.32.3",
]
all = ["mdvtools[app]"]           # alias

[project.scripts]
mdvtools = "mdvtools.cli:cli"

# dev/docs deps → PEP 735 [dependency-groups] (uv-native), NOT shipped to users
```

Changes that follow:
- **Build backend → `uv_build`; lockfile/installer → `uv`** (re-land PR #468). Poetry dropped. (§5, ADR-0002)
- **`mdvtools_lite_build/` is deleted** — no more `force-include`, second version, or second backend.
- **All modules ship** in the wheel (Path 1) — no physical exclusion. (ADR-0001)

---

## 4. Import sweep — what makes one-package-slim-by-default safe

**Headline findings (verified against the code)**
- `cli.py`, `mdvproject.py`, `conversions.py`, `serverlite.py`, `server_utils.py` import **zero** optional deps at top level. → `import mdvtools` + `MDVProject().serve()` survive a slim install untouched.
- **No core module imports the optional cluster at top level.** Arrows point only inward.
- **`spatial` is in core** (spatialdata required) → no spatial guards needed.
- **`auth` + `dbutils` + `llm` are a mutual import cycle** → modelled as the single `app` extra:

```
auth    → dbutils (project_auth → dbservice)   ┐
dbutils → auth    (mdv_server_app → authutils) │  one cycle → one `app` extra
dbutils → llm     (server_options → chat ext)  │  (cannot be installed independently)
llm     → auth    (chat_server_extension)      ┘
```

**6 leaf modules need an import guard** → all point to `mdvtools[app]`, via a shared `mdvtools/_optional.py` helper (`require_extra("app", <probe>)`):

| Package | File | Heavy top-level import | Probe |
|---|---|---|---|
| llm | `langchain_mdv.py` | langchain*, faiss | `langchain_openai` |
| llm | `chatlog.py` | langchain_core.* | `langchain_core` |
| llm | `github_utils.py` | requests, nbformat, dotenv | `nbformat` |
| dbutils | `dbmodels.py` | flask_sqlalchemy | `flask_sqlalchemy` |
| dbutils | `mdv_server_app.py` | sqlalchemy, psycogreen.gevent | `sqlalchemy` |
| auth | `auth0_provider.py` | authlib, jose, auth0.* | `authlib` |

(`dbservice.py`, `project_manager_extension.py`, `server_options.py` are heavy only *transitively*; guarding the leaves propagates a friendly error up the chain. **`llm/code_execution.py` needs no guard** — its only `matplotlib` reference is inside the `plt_code` *string* run in a subprocess; its own top-level imports are all stdlib.)

**Dependency-classification fixes (verified against runtime imports — as implemented)**
1. **`ruff` is in core but is a *runtime* dep** → move to **`app`** (and keep in **`dev`**). The chat feature shells out to `ruff check` (`llm/code_manipulation.py` `subprocess.run(['ruff', …])`), so it is needed at runtime; it is *also* a dev tool (`make lint`/`make format`). *(Earlier drafts said "core → dev, imported nowhere" — wrong: the CLI subprocess use was missed because nothing does `import ruff`.)*
2. **`generate-tiff-offsets` is in the dev group but imported by `spatial/jp2k.py`** → move to **core** (spatial is core).
3. **`matplotlib` is in the dev group but imported at runtime by `llm/langchain_mdv.py`** → move to **`app`**.
4. **`regex`, `python-dotenv`, `nbformat`, `tabulate`** (core) and **`gspread`, `oauth2client`** (the PEP 735 `llm` group) are llm-only → move to **`app`**. `langchain-community` (imported directly by `langchain_mdv.py`, previously only transitive) is declared explicitly in **`app`**.
5. dev-only confirmed clean: `leidenalg`, `umap-learn`, `seaborn`, `statsmodels`, `pillow`, `ome-types`, `spatialdata-plot`, `dummy-spatialdata`, `jupyterlab`, `pyright`, `pytest`, `setuptools`, `packaging`.

> **Note on `gspread`/`oauth2client`:** the old Dockerfile synced `--group dev --group backend --group auth` (never `llm`), so these were never installed in the image. Folding them into `app` means the full image now includes them — a (harmless) superset.

---

## 5. Build backend, versioning & uv

- **One version string** in `[project].version`. The two-version drift point disappears.
- **Backend/tooling = uv** (ADR-0002). uv already landed on `main` via **PR #472** (static PEP 621 metadata, `uv_build`, `uv.lock`, pinned uv, uv Dockerfile/CI). This work builds on top of it — we only add the `app` extra + guards. `uv_build` suffices because Path 1 needs no `force-include`.
- **Re-lock required:** moving deps between core / `app` / `dev` invalidates `uv.lock`; run `uv lock` (in the devcontainer at the pinned uv version) and commit the result. Verified: the re-lock changed **no versions and no packages** — pure re-classification.

### 5.1 Dockerfile impact

#472 already put the Dockerfile on uv. The only change here is the two `uv sync` lines: `--group dev --group backend --group auth` → `--group dev --extra app` (both the dependency-cache layer and the final install layer). `--extra app` pulls the whole server/chat/auth bundle. Everything else (in-project `.venv`, two-stage caching, `uv run -- gunicorn …`) is untouched.

---

## 6. Migration steps

1. ~~Recover why #471 reverted #468~~ — **moot:** uv re-landed on `main` via **PR #472**; we build on it.
2. ~~Re-land the uv work~~ — **done by #472** (static PEP 621, `uv_build`, `uv.lock`, uv Dockerfile/CI). Merged into this branch.
3. **Audit slim core** — confirmed `import mdvtools` + `MDVProject` + `.serve()` (serverlite/cli/conversions) import no optional deps at top level. ✅
4. **Reconcile disagreeing pins** (`polars`, `spatialdata`, `pandas`) — resolved by adopting #472's core pins as-is. ✅
5. **Add the single `app` extra** + `all` alias; apply the §4 classification fixes (`ruff`→app+dev, `generate-tiff-offsets`→core, `matplotlib`/`regex`/`dotenv`/`nbformat`/`tabulate`/`gspread`/`oauth2client`/`langchain-community`→app); delete the `backend`/`llm`/`auth` groups. ✅
6. **Add `mdvtools/_optional.py`** + guards on the 6 leaf modules (§4) + `tests/test_optional.py`. ✅
7. **Update the Dockerfile** — two `uv sync` lines → `--group dev --extra app` (§5.1). ✅
8. **Re-lock** — `uv lock` regenerated (no version/package drift); `uv sync --frozen --group dev --extra app` clean; `make test` green (186 passed). ✅
9. **Still to do:** confirm a *slim* env (`uv sync --no-default-groups`, or `pip install .` in a clean venv) imports + runs the CLI, and that a chat/server/auth call there raises the friendly `mdvtools[app]` error.
10. **Delete `mdvtools_lite_build/`** — done; its salvageable content was migrated: the PyPI quickstart → `python/PYPI_README.md` (now the package's PyPI long description, `readme = "PYPI_README.md"`), and the publish workflow → `docs/RELEASING.md` (rewritten for `uv build`/`uv publish` instead of Hatchling/twine). ✅
11. **Publish to TestPyPI first**, then PyPI — see `docs/RELEASING.md`. Release note: slim install is now lighter; full-app users need `mdvtools[app]`; spatial unchanged.

---

## 7. Open questions for standup

1. ~~Why was uv #468 reverted (#471)?~~ **Moot** — uv re-landed on `main` via #472 and is passing CI.
2. **`required-version` exact pin** — #472 set `[tool.uv] required-version = "0.11.14"` (exact `==`), which fails on any uv drift (hit during this work). Consider relaxing to a range (e.g. `>=0.11.14,<0.12`).
3. **PyPI ownership** — add Marouf (`mrfshk`) as collaborator on PyPI + TestPyPI `mdvtools` (Martin to action).
4. **Existing-user impact** — slim install becomes *lighter*; full-app users must add `[app]`. spatial is unchanged (stays in core). The full Docker image now also pulls `gspread`/`oauth2client` (previously absent). Note in release notes / version bump.

---

## 8. Appendix — approaches considered and not chosen

- **Keep Poetry** (the earlier recommendation). Poetry 2.x supports PEP 621 extras and would be a ~2-line Dockerfile change, but it leaves us on Poetry against the team's direction and re-converts the toml later anyway. Rejected — see ADR-0002.
- **Spatial as an extra** (Martin's preference). Today's PyPI wheel already ships spatial + requires spatialdata; demoting it breaks existing `pip install mdvtools` users for a maintainability goal that doesn't require the break. Rejected — spatial stays in core. (ADR-0001)
- **Separate named extras (spatial/chat/server/auth).** auth+dbutils+llm form an import cycle, so three "independent" extras that each pull the other two is false granularity. Collapsed to one `app` extra. (ADR-0001 / §4)
- **Two PyPI packages (pydantic-ai umbrella) / separate repo (Martin's "MDVApp").** Both *add* packages/versions to manage — the opposite of the goal. The 576 KB figure removes the footprint rationale. Rejected. (ADR-0001)

## 9. References

- pydantic-ai, Hugging Face transformers, pandas, scanpy — single-package, ship-all-modules + guarded imports.
- PEP 621 (`[project]`), PEP 735 (`[dependency-groups]`); uv docs (`uv_build`, `uv sync --all-extras`).
- ADR-0001 (single package + extras over module separation), ADR-0002 (adopt uv), `CONTEXT.md` (glossary).
