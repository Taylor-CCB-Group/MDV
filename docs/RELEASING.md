# Releasing `mdvtools` to PyPI

`mdvtools` is one package with one version, set in `python/pyproject.toml`.
`pip install mdvtools` installs the slim core, and `pip install "mdvtools[app]"` adds the full
application. See `MDVTOOLS_PACKAGING_PROPOSAL.md`, `docs/adr/0001` and `docs/adr/0002`.

GitHub Actions releases `mdvtools` every 4 months. A person reviews the changelog and approves
the PyPI upload, and the workflows do the rest. Every release goes to TestPyPI first, and the
workflow tests the TestPyPI upload before anything reaches PyPI. The manual steps at the end
are a fallback for when the workflows cannot run.

## Write changelog entries as you go

Each pull request that changes what users see adds a line under `## [Unreleased]` in
`CHANGELOG.md`, in the section that fits: `Added`, `Changed` or `Fixed`. Write it for someone
who uses MDV, and end it with the pull request number, such as `(#575)`. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

The release pull request turns `[Unreleased]` into the new version's section, so the entries
written during the 4 months become the release notes.

## How a release happens

1. At 9am UTC on 1 February, 1 June and 1 October, the Prepare release workflow
   (`.github/workflows/release-prepare.yml`) opens a pull request named "Release X.Y.Z". It
   bumps the minor version with `uv version --bump minor` and dates the changelog. If nothing
   has merged since the last release, it does nothing. If `[Unreleased]` is empty, it lists
   the titles of the pull requests merged since the last release under `### To rewrite` in
   the new section.
2. While the release pull request is open, each push to `main` runs the same workflow's
   `sync` job. It merges `main` into the release branch and adds the titles of the newly
   merged pull requests to the `### To rewrite` list. If the merge conflicts, it comments on
   the release pull request with the titles instead, and someone merges `main` into the
   branch by hand.
3. A maintainer reviews the release pull request:
   - Approve its waiting CI runs. GitHub holds checks on a pull request that a workflow
     opened until someone with write access approves them.
   - Edit the new `CHANGELOG.md` section so it describes the release for users. Rewrite the
     titles under `### To rewrite` into `Added`, `Changed` or `Fixed` entries, then delete
     the heading. The Release workflow stops while the heading is there.
   - Change the version if the release needs a patch or major bump instead. Follow
     [Semantic Versioning](https://semver.org/spec/v2.0.0.html): major for changes that
     break existing use, minor for new features, patch for fixes only.
4. Merging the pull request changes the version in `python/pyproject.toml` on `main`, which
   starts the Release workflow (`.github/workflows/release.yml`). It:
   - builds the frontend and the sdist and wheel
   - publishes them to TestPyPI
   - downloads the wheel back from TestPyPI, checks it matches the build, and smoke-tests it
     in a slim install and an `[app]` install (`.github/scripts/smoke_test.py`)
   - waits for a maintainer to approve the `pypi` environment
   - publishes to PyPI
   - tags `vX.Y.Z` and creates the GitHub release, with the changelog section as its notes
5. Approve the PyPI step from the workflow run once the TestPyPI checks have passed. Install
   from TestPyPI yourself first if the release needs more testing than the smoke test does.

To release outside the schedule, run the Prepare release workflow from the Actions tab and
choose the bump. You can also open a pull request yourself that bumps the version with
`uv version --bump` (from `python/`) and dates the changelog with
`python .github/scripts/release.py changelog --previous <old version>`. The Release workflow
runs whenever a version without a tag reaches `main`, however the bump got there.

## One-off setup

The workflows publish with PyPI trusted publishing, so no API tokens are stored in GitHub.
This setup needs an admin of the GitHub repository and an owner of `mdvtools` on PyPI and
TestPyPI.

1. On pypi.org, open "Your projects", choose Manage for `mdvtools`, then Publishing. Add a
   GitHub Actions publisher with owner `Taylor-CCB-Group`, repository `MDV`, workflow
   `release.yml` and environment `pypi`.
2. Do the same on test.pypi.org, with environment `testpypi`. TestPyPI has separate accounts
   from PyPI.
3. In the GitHub repository, go to Settings, then Environments. Create `testpypi` and `pypi`,
   and add required reviewers to `pypi`. Those reviewers approve each PyPI upload.
4. In Settings, then Actions, then General, turn on "Allow GitHub Actions to create and
   approve pull requests", so the Prepare release workflow can open its pull request.

## Release by hand

Use these steps when the workflows cannot run, such as when GitHub Actions is down. Run every command from the repository root.

### Prerequisites

- accounts on TestPyPI and PyPI, with access to the `mdvtools` project on each
- an API token for each index, exported as `UV_PUBLISH_TOKEN` before you publish
- `uv` at the version pinned in `[tool.uv] required-version`
- `pnpm` for the frontend build

### 1. Bump the version and date the changelog

On a new branch:

```bash
previous=$(python .github/scripts/release.py current)
(cd python && uv version --bump minor --no-sync)    # or patch or major
python .github/scripts/release.py changelog --previous "$previous"
```

Commit these changes and merge them into `main` through a pull request, then check out the
merged `main`. You build and tag from that commit, so the tag matches the published package:

```bash
git switch main && git pull
```

If GitHub Actions is working again when the pull request merges, the Release workflow
publishes the release itself, and you can stop here.

### 2. Build the frontend

The wheel ships the compiled JavaScript and CSS from `python/mdvtools/static`, which git
ignores. Build it before you package:

```bash
pnpm install
pnpm run build-flask-vite        # writes to python/mdvtools/static
```

### 3. Build the sdist and wheel

```bash
(cd python && uv build)          # writes the sdist and wheel to python/dist/
```

Check that the wheel includes `mdvtools/static` and leaves out notebooks and `*.map` files:

```bash
python -m zipfile -l python/dist/mdvtools-*.whl | grep -E "static/|\.ipynb$" | head
```

### 4. Publish to TestPyPI and test it

```bash
# UV_PUBLISH_TOKEN holds your TestPyPI token
uv publish --publish-url https://test.pypi.org/legacy/ python/dist/*
```

Install the wheel in clean environments and run the smoke test against both installs:

```bash
pip download --no-deps --only-binary :all: --index-url https://test.pypi.org/simple/ \
  -d /tmp/testpypi "mdvtools==X.Y.Z"
uv venv --python 3.12 /tmp/mdv-slim && uv pip install --python /tmp/mdv-slim /tmp/testpypi/*.whl
/tmp/mdv-slim/bin/python .github/scripts/smoke_test.py X.Y.Z
/tmp/mdv-slim/bin/mdvtools --version   # should report X.Y.Z and the release commit
uv venv --python 3.12 /tmp/mdv-app && uv pip install --python /tmp/mdv-app "$(ls /tmp/testpypi/*.whl)[app]"
/tmp/mdv-app/bin/python .github/scripts/smoke_test.py X.Y.Z --app
```

Installing the downloaded wheel, rather than installing from the TestPyPI index, makes sure
every dependency comes from PyPI.

### 5. Publish to PyPI

```bash
# UV_PUBLISH_TOKEN now holds your PyPI token
uv publish python/dist/*
```

### 6. Tag the release and publish the notes

```bash
python .github/scripts/release.py notes X.Y.Z > /tmp/notes.md
gh release create vX.Y.Z python/dist/* --target "$(git rev-parse HEAD)" \
  --title "MDV X.Y.Z" --notes-file /tmp/notes.md
```

`gh release create` creates the tag on the commit you built from in step 1. A tagged version
stops the Release workflow from publishing it again.

## Notes

- `[tool.uv] required-version` pins uv to one exact version, and any other version fails
  `uv build`, `uv lock` and `uv sync`. Install the pinned version, or run it with
  `uvx uv@<version>`.
- The PyPI page shows `python/PYPI_README.md`, not the developer `README.md`. Update its
  quickstart when the API changes. `CHANGELOG.md` appears on GitHub only.
- `uv build` and `uv publish` do not read `uv.lock`. The lock matters for `uv sync` in
  development, Docker and CI.
- After each `sync` push to the release branch, the new CI runs on the release pull request
  wait for someone with write access to approve them. GitHub holds every pull request run
  that a `GITHUB_TOKEN` push starts.
- PyPI and TestPyPI never accept the same file twice. If a release fails after its TestPyPI
  upload, bump to the next patch version rather than rebuilding the same one.
