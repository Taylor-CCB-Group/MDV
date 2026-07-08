import requests
from flask import Flask
from unittest.mock import patch

from mdvtools.ucsc_proxy_extension import UcscProxyServerExtension


class FakeUCSCResponse:
    def __init__(self, content: bytes, content_type: str, *, raise_for_status: bool = False):
        self._content = content
        self.headers = {"Content-Type": content_type}
        self._raise_for_status = raise_for_status
        self.closed = False

    def raise_for_status(self):
        if self._raise_for_status:
            raise requests.exceptions.RequestException("upstream error")

    def close(self):
        self.closed = True

    def iter_content(self, chunk_size=8192):
        for i in range(0, len(self._content), chunk_size):
            yield self._content[i : i + chunk_size]


def test_ucsc_proxy_allows_host_and_returns_image():
    app = Flask(__name__)
    app.config["MAX_UCSC_BYTES"] = 1024
    UcscProxyServerExtension().register_global_routes(app, app.config)

    content = b"fake-png-bytes"
    fake = FakeUCSCResponse(content, "image/png")

    with patch("mdvtools.ucsc_proxy_extension.requests.get", return_value=fake) as get:
        resp = app.test_client().get(
            "/ucsc_proxy",
            query_string={"ucscHost": "genome.ucsc.edu", "db": "hg38", "position": "chr1:1-2"},
        )

    assert resp.status_code == 200
    assert resp.data == content
    assert get.call_count == 1


def test_ucsc_proxy_rejects_unknown_host():
    app = Flask(__name__)
    UcscProxyServerExtension().register_global_routes(app, app.config)

    with patch("mdvtools.ucsc_proxy_extension.requests.get") as get:
        resp = app.test_client().get(
            "/ucsc_proxy",
            query_string={"ucscHost": "evil.example.com", "db": "hg38"},
        )

    assert resp.status_code == 400
    assert b"Invalid host" in resp.data
    assert get.call_count == 0


def test_ucsc_proxy_rejects_non_image_content_type():
    app = Flask(__name__)
    UcscProxyServerExtension().register_global_routes(app, app.config)

    fake = FakeUCSCResponse(b"not-an-image", "text/plain")

    with patch("mdvtools.ucsc_proxy_extension.requests.get", return_value=fake):
        resp = app.test_client().get(
            "/ucsc_proxy",
            query_string={"ucscHost": "genome.ucsc.edu", "db": "hg38"},
        )

    assert resp.status_code == 415
    assert b"Unsupported media type" in resp.data


def test_ucsc_proxy_enforces_size_limit():
    app = Flask(__name__)
    app.config["MAX_UCSC_BYTES"] = 10
    UcscProxyServerExtension().register_global_routes(app, app.config)

    # total content (15 bytes) should exceed the 10-byte limit.
    fake = FakeUCSCResponse(b"0123456789abcde", "image/png")

    with patch("mdvtools.ucsc_proxy_extension.requests.get", return_value=fake):
        resp = app.test_client().get(
            "/ucsc_proxy",
            query_string={"ucscHost": "genome.ucsc.edu", "db": "hg38"},
        )

    assert resp.status_code == 502
    assert b"Upstream response too large" in resp.data


def test_ucsc_proxy_returns_502_on_upstream_failure():
    app = Flask(__name__)
    UcscProxyServerExtension().register_global_routes(app, app.config)

    with patch(
        "mdvtools.ucsc_proxy_extension.requests.get",
        side_effect=requests.exceptions.RequestException("boom"),
    ):
        resp = app.test_client().get(
            "/ucsc_proxy",
            query_string={"ucscHost": "genome.ucsc.edu", "db": "hg38"},
        )

    assert resp.status_code == 502
    assert b"Upstream fetch failed" in resp.data

