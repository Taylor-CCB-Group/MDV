"""A browser that visited one deploy must run the next deploy's code, not a cached copy.

Models what bit OOD users (kir-rescomp/ood-app-mdv#36): the app is redeployed under the same
URL, and the browser (or something between it and Flask) still holds files from before.
The real templates, `entry_assets` and `add_safe_headers` serve the pages; `Browser` is a
small HTTP cache (no-store, max-age, validators, heuristic freshness) that loads a page the
way a browser would: HTML, then the entry's static imports, then its lazy imports.

Not defended, and not tested: a proxy that caches the HTML itself. The page then keeps
naming the old entry whatever the build calls its files.
"""

import json
import os
import re
import time
from dataclasses import dataclass, field
from email.utils import parsedate_to_datetime
from urllib.parse import urljoin

import pytest
from flask import Flask, render_template

from mdvtools.entry_assets import register_entry_assets
from mdvtools.server_utils import add_safe_headers

TEMPLATES = os.path.join(os.path.dirname(os.path.dirname(__file__)), "templates")
HOST = "http://ood.example.org"
PAGE = f"{HOST}/project/1/"
MONTH = 30 * 24 * 3600


def write(path, text, mtime):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    os.utime(path, (mtime, mtime))


def make_build(folder, tag: str, *, hashed: bool, mtime: float | None = None):
    """A tiny stand-in for a frontend build whose entry names `tag` when it runs.

    Chunks are always hashed (as Vite does); the entry is only when `hashed`, which also
    writes the manifest. Without `hashed` it is the pre-#601 layout: js/mdv.js.
    """
    mtime = mtime if mtime is not None else time.time() - MONTH
    entry, css = (f"js/mdv-{tag}.js", f"assets/mdv-{tag}.css") if hashed else ("js/mdv.js", "assets/mdv.css")
    write(
        folder / entry,
        f'import "../assets/shared-{tag}.js"; const lazy = () => import("../assets/ChartManager-{tag}.js"); console.log("build {tag}");',
        mtime,
    )
    write(folder / css, f"/* {tag} */", mtime)
    write(folder / f"assets/shared-{tag}.js", f"/* {tag} */", mtime)
    write(folder / f"assets/ChartManager-{tag}.js", f"/* {tag} */", mtime)
    if hashed:
        manifest = {"src/modules/project_bootstrap.tsx": {"file": entry, "name": "mdv", "isEntry": True, "css": [css]}}
        write(folder / "manifest.json", json.dumps(manifest), mtime)
    return folder


class Proxy:
    """WSGI middleware standing in for a reverse proxy that may override Cache-Control."""

    def __init__(self, wsgi):
        self.wsgi = wsgi
        self.caches_assets = False

    def __call__(self, environ, start_response):
        def respond(status, headers, exc_info=None):
            if self.caches_assets and environ["PATH_INFO"].startswith("/flask/"):
                kept = [(k, v) for k, v in headers if k.lower() not in ("cache-control", "etag", "last-modified")]
                headers = [*kept, ("Cache-Control", "public, max-age=86400")]
            return start_response(status, headers, exc_info)

        return self.wsgi(environ, respond)


@dataclass
class Visit:
    ran: list[str] = field(default_factory=list)  # tags of the entries that executed
    failed: list[str] = field(default_factory=list)  # "404 /url" for anything that didn't load


@dataclass
class Stored:
    headers: dict
    body: str
    stored_at: float


class Browser:
    def __init__(self, client):
        self.client = client
        self.cache: dict[str, Stored] = {}
        self.now = time.time()

    def advance(self, seconds):
        self.now += seconds

    def _fresh(self, entry: Stored) -> bool:
        cc = entry.headers.get("Cache-Control", "")
        if "no-cache" in cc or "must-revalidate" in cc:
            return False
        age = self.now - entry.stored_at
        if m := re.search(r"max-age=(\d+)", cc):
            return age < int(m[1])
        if "Last-Modified" in entry.headers and "Date" in entry.headers:
            modified = parsedate_to_datetime(entry.headers["Last-Modified"])
            date = parsedate_to_datetime(entry.headers["Date"])
            return age < 0.1 * (date - modified).total_seconds()  # RFC 9111 heuristic
        return False

    def get(self, url: str) -> tuple[int, str]:
        path = url.removeprefix(HOST)
        entry = self.cache.get(url)
        if entry and self._fresh(entry):
            return 200, entry.body
        validators = {}
        if entry and "ETag" in entry.headers:
            validators["If-None-Match"] = entry.headers["ETag"]
        if entry and "Last-Modified" in entry.headers:
            validators["If-Modified-Since"] = entry.headers["Last-Modified"]
        response = self.client.get(path, headers=validators)
        if response.status_code == 304 and entry:
            self.cache[url] = Stored(dict(response.headers) | {k: v for k, v in entry.headers.items() if k not in response.headers}, entry.body, self.now)
            return 200, entry.body
        body = response.get_data(as_text=True)
        if response.status_code == 200 and "no-store" not in response.headers.get("Cache-Control", ""):
            self.cache[url] = Stored(dict(response.headers), body, self.now)
        return response.status_code, body

    def visit(self, url: str = PAGE) -> Visit:
        """Load the page, run its entry, then let the entry load its lazy chunks."""
        visit = Visit()
        status, html = self.get(url)
        assert status == 200
        for ref in re.findall(r'(?:href|src)="([^"]+)"', html):
            self._load(urljoin(url, ref), visit)
        return visit

    def _load(self, url: str, visit: Visit, seen: set[str] | None = None):
        seen = seen if seen is not None else set()
        if url in seen:
            return
        seen.add(url)
        status, body = self.get(url)
        if status != 200:
            visit.failed.append(f"{status} {url.removeprefix(HOST)}")
            return
        for ref in re.findall(r'import\s*\(?\s*"([^"]+)"', body):
            self._load(urljoin(url, ref), visit, seen)
        if m := re.search(r'console\.log\("build (\w+)"\)', body):
            visit.ran.append(m[1])


@pytest.fixture
def site(tmp_path):
    """A Flask app whose static folder is swapped for each 'deploy', behind a Proxy."""
    app = Flask(__name__, template_folder=TEMPLATES, static_folder=str(tmp_path / "empty"), static_url_path="/flask")
    app.url_map.strict_slashes = False
    app.after_request(add_safe_headers)
    register_entry_assets(app)

    @app.route("/project/<project_id>/")
    def project(project_id):
        return render_template("page.html")

    proxy = Proxy(app.wsgi_app)
    app.wsgi_app = proxy

    class Site:
        browser = Browser(app.test_client())

        def deploy(self, folder):
            app.static_folder = str(folder)

        caches_assets = property(lambda self: proxy.caches_assets, lambda self, v: setattr(proxy, "caches_assets", v))

    return Site()


@pytest.mark.parametrize("proxy_caches_assets", [False, True])
@pytest.mark.parametrize("before_hashed", [True, False], ids=["hashed->hashed", "pre-#601->hashed"])
def test_next_deploy_runs_new_code(site, tmp_path, before_hashed, proxy_caches_assets):
    site.caches_assets = proxy_caches_assets
    site.deploy(make_build(tmp_path / "a", "AAAA1111", hashed=before_hashed))
    first = site.browser.visit()
    assert first.ran == ["AAAA1111"] and not first.failed

    site.browser.advance(3600)
    site.deploy(make_build(tmp_path / "b", "BBBB2222", hashed=True))
    second = site.browser.visit()

    assert second.ran == ["BBBB2222"]
    assert second.failed == []


def test_unhashed_entries_are_revalidated_by_an_honest_server(site, tmp_path):
    # the fallback for builds without a manifest still has to notice a new deploy
    site.deploy(make_build(tmp_path / "a", "AAAA1111", hashed=False, mtime=time.time() - MONTH))
    site.browser.visit()

    site.browser.advance(3600)
    site.deploy(make_build(tmp_path / "b", "BBBB2222", hashed=False, mtime=time.time() - 60))

    assert site.browser.visit().ran == ["BBBB2222"]


def test_unhashed_entries_behind_a_caching_proxy_go_stale(site, tmp_path):
    # what the manifest exists to prevent, and why a hashed build without one is an error
    site.caches_assets = True
    site.deploy(make_build(tmp_path / "a", "AAAA1111", hashed=False))
    site.browser.visit()

    site.browser.advance(3600)
    site.deploy(make_build(tmp_path / "b", "BBBB2222", hashed=False))

    assert site.browser.visit().ran == ["AAAA1111"]
