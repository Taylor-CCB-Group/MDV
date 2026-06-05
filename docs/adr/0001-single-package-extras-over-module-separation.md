# Single `mdvtools` package with optional-dependency extras, not physical module separation

**Status:** accepted

## Context & decision

`mdvtools` historically shipped from **two** diverging `pyproject.toml` files: the full
Poetry build in `python/` (drives the Docker image, never published) and a slim Hatchling
build in `mdvtools_lite_build/` that used `force-include` to publish a hand-picked *subset*
of modules to PyPI. The two drifted (disagreeing pins, two version strings, a stale
include list). The goal — set by management — is to **manage one `pyproject.toml`**, not
to maintain two.

We will collapse to a **single `python/pyproject.toml`**: required dependencies = the slim
core, and heavy optional features become `[project.optional-dependencies]` extras guarded
by friendly import errors. `pip install mdvtools` is slim by default; the Docker image
installs all extras for the full experience. **`mdvtools_lite_build/` is deleted.** All
modules ship in the wheel (Path 1); we do **not** physically exclude `dbutils`/`auth`/`llm`
from the published wheel.

## Why physical module separation was rejected

The lite build's `force-include` trims `.py` *modules*, but the install footprint is
dominated by *dependencies*. Martin (who built the lite wheel) measured this directly:
cherry-picking modules with Hatchling **still produced a >1.4 GB install**, because the
weight is langchain/faiss/scanpy/spatialdata et al., not source files. Excluding all of
`dbutils`+`auth`+`llm`+`spatial` source saves ~576 KB out of ~1.2–1.4 GB (~0.05%) — and the
59 MB built-JS bundle ships in every wheel regardless. So module-level file selection
operates on the axis that doesn't matter; extras control the axis that does, and are
backend-agnostic.

This also matches ecosystem norms. The two established patterns are **(A)** one package
that ships all modules and gates features with optional-deps + guarded imports
(pandas, HuggingFace transformers, scanpy, scikit-learn), and **(B)** genuinely separate
packages each with its own source tree (langchain-core / -community / -openai). Publishing
a hand-picked subset of one package's modules via the build backend is neither — it is a
home-grown third option that no major project maintains, and it is brittle (stale include
list, published wheel ≠ editable/source install). We adopt Pattern A.

## Considered and rejected

- **Path 2 — physically exclude app modules from the slim wheel.** Buys ~0.05% size for a
  permanent maintenance burden and a published-wheel-≠-source footgun. Also awkward under
  `uv_build`. Rejected.
- **Separate repository (Martin's "MDVApp").** Largest restructuring cost (moves import
  paths, own release cadence); justified only by hard IP/separation needs we don't have.
  The 576 KB figure removes the footprint rationale. Rejected.

## Consequences

- The build backend no longer needs surgical module selection; `uv_build`'s basic
  include/exclude (for notebooks, scratch JSON, test data, `*.map`) suffices.
- Leaf modules with heavy top-level imports need guarded imports so a slim install fails
  with "install `mdvtools[<extra>]`" rather than an `ImportError`.
- `spatialdata` stays in the **slim core** (not an extra), preserving current PyPI
  behaviour — a deliberate deviation from Martin's "spatial as an extra" preference, since
  the maintainability goal does not require breaking existing `pip install mdvtools` users.
