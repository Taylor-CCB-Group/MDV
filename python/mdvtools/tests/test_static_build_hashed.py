"""Every script and stylesheet in a built frontend has its content hash in its name.

A file served under a fixed name can change between deploys, and a copy cached before the
change then runs against code that no longer matches it (kir-rescomp/ood-app-mdv#36).
Checks the Flask build in `MDV_STATIC_DIR` (default: mdvtools/static, from
`pnpm run build-flask-vite`); skipped when that isn't a hashed build.
"""

import os
import re
from pathlib import Path

import pytest

STATIC = Path(os.environ.get("MDV_STATIC_DIR", Path(__file__).parents[1] / "static"))
HASHED = re.compile(r"-[A-Za-z0-9_-]{8}\.(js|css)$")
BUILT_DIRS = ("js", "assets", "img")


@pytest.mark.skipif(not (STATIC / "manifest.json").exists(), reason=f"no hashed frontend build in {STATIC}")
def test_scripts_and_styles_are_content_hashed():
    unhashed = [
        str(path.relative_to(STATIC))
        for folder in BUILT_DIRS
        for path in (STATIC / folder).rglob("*")
        if path.suffix in (".js", ".css") and not HASHED.search(path.name)
    ]

    assert unhashed == []
