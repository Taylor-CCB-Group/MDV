from flask import Flask, Response

from mdvtools.server_utils import add_safe_headers


def create_test_app():
    app = Flask(__name__)
    app.after_request(add_safe_headers)

    @app.route("/")
    @app.route("/project/<project_id>/")
    def page(project_id=None):
        return Response("<!doctype html><div id='holder'></div>", mimetype="text/html")

    @app.route("/static/js/<path:name>")
    def js(name):
        response = Response("console.log('ok')", mimetype="text/javascript")
        response.headers["Cache-Control"] = "public, max-age=31536000"
        return response

    @app.route("/static/assets/<path:name>")
    def asset(name):
        response = Response("body {}", mimetype="text/css")
        response.headers["Cache-Control"] = "public, max-age=31536000"
        return response

    @app.route("/state.json")
    def state_json():
        return Response("{}", mimetype="application/json")

    @app.route("/images/<path:name>")
    def image(name):
        response = Response(b"image", mimetype="image/png")
        response.headers["Cache-Control"] = "public, max-age=31536000"
        return response

    return app


def test_html_shells_are_not_stored():
    client = create_test_app().test_client()

    assert client.get("/").headers["Cache-Control"] == "no-store"
    assert client.get("/project/abc/").headers["Cache-Control"] == "no-store"


def test_stable_vite_entry_assets_revalidate():
    client = create_test_app().test_client()

    for path in [
        "/static/js/mdv.js",
        "/static/js/catalog.js",
        "/static/js/login.js",
        "/static/assets/mdv.css",
        "/static/assets/catalog.css",
    ]:
        assert (
            client.get(path).headers["Cache-Control"]
            == "no-cache, max-age=0, must-revalidate"
        )


def test_project_data_and_explicit_cache_headers_are_not_overridden():
    client = create_test_app().test_client()

    assert "Cache-Control" not in client.get("/state.json").headers
    assert (
        client.get("/static/assets/ChartManager-Abc123.js").headers["Cache-Control"]
        == "public, max-age=31536000"
    )
    assert (
        client.get("/images/example.png").headers["Cache-Control"]
        == "public, max-age=31536000"
    )
