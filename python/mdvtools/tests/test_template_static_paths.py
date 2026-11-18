"""The HTML shells must reach their static assets behind a prefix-stripping proxy.

Some proxies, Verily Workbench among them, serve MDV under a path such as
`/app/<uuid>/proxy/5055/` and strip that prefix before the request reaches Flask.
Flask then sees `/project/7/` and `url_for('static')` builds `/flask/js/mdv.js`,
which the browser resolves against the proxy host's root and misses MDV entirely.
"""

import os
import re
from urllib.parse import urljoin

import pytest
from flask import Flask, render_template

from mdvtools.entry_assets import register_entry_assets

TEMPLATES = os.path.join(os.path.dirname(os.path.dirname(__file__)), "templates")
PROXY_ROOT = "https://workbench.example.org/app/abc/proxy/5055/"


def create_test_app(static_folder):
    # Same static URL as create_flask_app in dbutils/mdv_server_app.py (/flask/..., not
    # Flask's default /static/...). The folder is an empty one from the test, not the real
    # build in /app/dist/flask: a hashed build there (with manifest.json) would change the
    # file names, and these tests are about paths, not which build is installed.
    app = Flask(
        __name__,
        template_folder=TEMPLATES,
        static_folder=str(static_folder),
        static_url_path="/flask",
    )
    app.url_map.strict_slashes = False
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

    return app


def asset_urls(html):
    return re.findall(r'(?:href|src)="([^"]+)"', html)


@pytest.mark.parametrize(
    "path, expected",
    [
        ("", ["flask/assets/catalog.css", "flask/js/catalog.js"]),
        ("login", ["flask/js/login.js"]),
        ("project/7/", ["flask/assets/mdv.css", "flask/js/mdv.js"]),
        ("project/7", ["flask/assets/mdv.css", "flask/js/mdv.js"]),
    ],
)
def test_assets_resolve_under_the_proxy_prefix(path, expected, tmp_path):
    html = create_test_app(tmp_path).test_client().get("/" + path).get_data(as_text=True)

    page_url = PROXY_ROOT + path
    resolved = [urljoin(page_url, url) for url in asset_urls(html)]

    assert resolved == [PROXY_ROOT + url for url in expected]


def test_assets_match_url_for_without_a_proxy(tmp_path):
    client = create_test_app(tmp_path).test_client()

    html = client.get("/project/7/", base_url="http://localhost/mdv/").get_data(as_text=True)
    resolved = [urljoin("http://localhost/mdv/project/7/", url) for url in asset_urls(html)]

    assert resolved == [
        "http://localhost/mdv/flask/assets/mdv.css",
        "http://localhost/mdv/flask/js/mdv.js",
    ]
