"""Find the built files for each frontend entry point (mdv, catalog, login).

The Flask build names its entries with a content hash (`js/mdv-<hash>.js`), so a cached
copy can't outlive a deploy and go on importing chunks that no longer exist. Vite records
the hashed names in `manifest.json` at the root of the static folder. Builds without a
manifest still use the old unhashed names.
"""

import json
import os
from dataclasses import dataclass
from functools import lru_cache

from flask import Flask

MANIFEST = "manifest.json"


@dataclass(frozen=True)
class EntryAssets:
    js: str
    css: tuple[str, ...]


UNHASHED = {
    "mdv": EntryAssets("js/mdv.js", ("assets/mdv.css",)),
    "catalog": EntryAssets("js/catalog.js", ("assets/catalog.css",)),
    "login": EntryAssets("js/login.js", ()),
}


def entry_assets(static_folder: str, name: str) -> EntryAssets:
    """Paths relative to the static folder for entry `name`."""
    path = os.path.join(static_folder, MANIFEST)
    try:
        # keyed by mtime so a rebuild under a running server is picked up
        mtime = os.stat(path).st_mtime_ns
    except FileNotFoundError:
        return UNHASHED[name]
    entries = _read_manifest(path, mtime)
    if name not in entries:
        raise LookupError(f"No entry '{name}' in {path}")
    return entries[name]


@lru_cache(maxsize=4)
def _read_manifest(path: str, mtime: int) -> dict[str, EntryAssets]:
    with open(path) as f:
        manifest = json.load(f)
    return {
        chunk["name"]: EntryAssets(chunk["file"], _css(manifest, key))
        for key, chunk in manifest.items()
        if chunk.get("isEntry") and "name" in chunk
    }


def _css(manifest: dict, key: str) -> tuple[str, ...]:
    # Vite injects CSS for dynamic imports itself; CSS of the entry and anything it imports
    # statically has to be linked from the page.
    css: list[str] = []
    seen: set[str] = set()

    def visit(k: str):
        if k in seen:
            return
        seen.add(k)
        chunk = manifest[k]
        for imported in chunk.get("imports", []):
            visit(imported)
        for file in chunk.get("css", []):
            if file not in css:
                css.append(file)

    visit(key)
    return tuple(css)


def register_entry_assets(app: Flask):
    """Expose `entry_assets(name)` to this app's templates."""

    def template_entry_assets(name: str) -> EntryAssets:
        if app.static_folder is None:
            return UNHASHED[name]
        return entry_assets(app.static_folder, name)

    app.add_template_global(template_entry_assets, "entry_assets")
