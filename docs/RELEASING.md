# Releasing `mdvtools` to PyPI

How to build and publish the `mdvtools` package. One package, one version, one
`pyproject.toml` (`python/pyproject.toml`) — `pip install mdvtools` is the slim core and
`pip install "mdvtools[app]"` adds the full app (see `MDVTOOLS_PACKAGING_PROPOSAL.md`,
`docs/adr/0001`, `docs/adr/0002`).

> Replaces the old `mdvtools_lite_build/` flow (Hatchling + `twine`). We now build and
> publish with **uv**.

## Prerequisites

- Accounts on **TestPyPI** and **PyPI**, and **collaborator** access on the `mdvtools`
  project on each (ask Martin Sergeant to add you).
- An **API token** for each index (account → API tokens). Export it when publishing:
  ```bash
  export UV_PUBLISH_TOKEN="pypi-…"     # the token for the index you're publishing to
  ```
- `uv` on your PATH at the version pinned in `[tool.uv] required-version` (CI/Docker use it).
- `pnpm` for the frontend build.

## Steps

### 1. Bump the version

Single source of truth — edit `[project].version` in `python/pyproject.toml`:

```toml
[project]
version = "1.2.4"     # bump per semver
```

Note for the changelog: the slim install is intentionally lighter than older releases —
full-app users must now install `mdvtools[app]`. Spatial support stays in the slim core.

### 2. Build the frontend assets (required)

The wheel ships the compiled JS/CSS from `python/mdvtools/static`. Build it **before**
packaging so the release contains current frontend code:

```bash
pnpm install
pnpm run build-flask-vite        # outputs into python/mdvtools/static
```

### 3. Build the distributions with uv

```bash
cd python
uv build                         # writes sdist + wheel to python/dist/
```

Sanity-check the wheel contents (should include `mdvtools/static`, exclude notebooks /
scratch JSON / `*.map`):

```bash
python -m zipfile -l dist/mdvtools-*.whl | grep -E "static/|\.ipynb$" | head
```

### 4. Publish to TestPyPI first

```bash
# UV_PUBLISH_TOKEN must hold your *TestPyPI* token
uv publish --publish-url https://test.pypi.org/legacy/ dist/*
```

Then verify in a clean environment — both the slim and full installs, plus a smoke test:

```bash
python -m venv /tmp/mdv-slim && . /tmp/mdv-slim/bin/activate
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ mdvtools
mdvtools --help
# a full-app feature here should raise a friendly: pip install "mdvtools[app]"
deactivate

python -m venv /tmp/mdv-app && . /tmp/mdv-app/bin/activate
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ "mdvtools[app]"
deactivate
```

(The `--extra-index-url` lets TestPyPI installs resolve real dependencies from PyPI.)

### 5. Publish to PyPI

Once TestPyPI looks good:

```bash
# UV_PUBLISH_TOKEN must now hold your *PyPI* token
cd python
uv publish dist/*
```

### 6. Tag the release

```bash
git tag -a v1.2.4 -m "mdvtools 1.2.4" && git push --tags
```

## Notes & gotchas

- **`required-version` is an exact pin** (`[tool.uv] required-version = "0.11.14"`). Any other
  uv version fails `uv build`/`uv lock`/`uv sync`. If you hit it, match the version
  (`pip install 'uv==<pinned>'`) and run `hash -r` (bash caches the `uv` path). Relaxing this
  to a range is a tracked follow-up.
- **The PyPI page renders `python/PYPI_README.md`** (`readme = "PYPI_README.md"`), not the
  dev-facing `README.md`. Update the quickstart there if the API changes.
- `uv build`/`uv publish` do **not** need `uv.lock`; the lock matters for `uv sync`
  (dev/Docker/CI), not for the published wheel's dependency metadata.
- **Fallback with twine** if needed: `uv build` then
  `twine upload --repository testpypi dist/*` / `twine upload dist/*`.
