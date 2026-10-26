"""The HTML shells link the hashed entry files that the Flask build records in manifest.json."""

import json
import os
import re

from flask import Flask, render_template

from mdvtools.entry_assets import entry_assets, register_entry_assets

TEMPLATES = os.path.join(os.path.dirname(os.path.dirname(__file__)), "templates")

# Trimmed from a real `build=desktop_pt` manifest.
MANIFEST = {
    "src/modules/project_bootstrap.tsx": {
        "file": "js/mdv-B4coGEvZ.js",
        "name": "mdv",
        "src": "src/modules/project_bootstrap.tsx",
        "isEntry": True,
        "imports": ["_shared-Xy12.js"],
        "dynamicImports": ["_ChartManager-BndN7OT2.js"],
        "css": ["assets/mdv-CVOPvuaJ.css"],
    },
    "src/catalog/catalog_index.tsx": {
        "file": "js/catalog-CjSaRnjV.js",
        "name": "catalog",
        "src": "src/catalog/catalog_index.tsx",
        "isEntry": True,
        "imports": ["_shared-Xy12.js"],
        "css": ["assets/catalog-C040d4dB.css"],
    },
    "src/login/login_index.tsx": {
        "file": "js/login-EahKsMqf.js",
        "name": "login",
        "src": "src/login/login_index.tsx",
        "isEntry": True,
    },
    "_shared-Xy12.js": {
        "file": "assets/shared-Xy12.js",
        "name": "shared",
        "css": ["img/shared-Ab34.css"],
    },
    "_ChartManager-BndN7OT2.js": {
        "file": "assets/ChartManager-BndN7OT2.js",
        "name": "ChartManager",
        "isDynamicEntry": True,
        "css": ["img/ChartManager-Cd56.css"],
    },
}


def write_manifest(folder, manifest=MANIFEST):
    with open(os.path.join(folder, "manifest.json"), "w") as f:
        json.dump(manifest, f)


def test_unhashed_names_without_a_manifest(tmp_path):
    assets = entry_assets(str(tmp_path), "mdv")

    assert assets.js == "js/mdv.js"
    assert assets.css == ("assets/mdv.css",)


def test_entry_css_includes_static_imports_but_not_dynamic_ones(tmp_path):
    write_manifest(tmp_path)

    assets = entry_assets(str(tmp_path), "mdv")

    assert assets.js == "js/mdv-B4coGEvZ.js"
    assert assets.css == ("img/shared-Ab34.css", "assets/mdv-CVOPvuaJ.css")
    assert entry_assets(str(tmp_path), "login").css == ()


def test_a_rebuild_is_picked_up(tmp_path):
    write_manifest(tmp_path)
    assert entry_assets(str(tmp_path), "login").js == "js/login-EahKsMqf.js"

    rebuilt = json.loads(json.dumps(MANIFEST))
    rebuilt["src/login/login_index.tsx"]["file"] = "js/login-New123.js"
    write_manifest(tmp_path, rebuilt)
    stat = os.stat(tmp_path / "manifest.json")
    os.utime(tmp_path / "manifest.json", ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))

    assert entry_assets(str(tmp_path), "login").js == "js/login-New123.js"


def test_templates_link_the_hashed_files(tmp_path):
    write_manifest(tmp_path)
    app = Flask(__name__, template_folder=TEMPLATES, static_folder=str(tmp_path), static_url_path="/flask")
    register_entry_assets(app)

    @app.route("/")
    def index():
        return render_template("index.html")

    @app.route("/login")
    def login():
        return render_template("login.html")

    @app.route("/project/<project_id>/")
    def project(project_id):
        return render_template("page.html")

    client = app.test_client()

    def asset_urls(path):
        return re.findall(r'(?:href|src)="([^"]+)"', client.get(path).get_data(as_text=True))

    assert asset_urls("/") == [
        "flask/img/shared-Ab34.css",
        "flask/assets/catalog-C040d4dB.css",
        "flask/js/catalog-CjSaRnjV.js",
    ]
    assert asset_urls("/login") == ["flask/js/login-EahKsMqf.js"]
    assert asset_urls("/project/7/") == [
        "../../flask/img/shared-Ab34.css",
        "../../flask/assets/mdv-CVOPvuaJ.css",
        "../../flask/js/mdv-B4coGEvZ.js",
    ]
