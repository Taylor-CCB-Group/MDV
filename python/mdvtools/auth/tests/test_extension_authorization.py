from flask import Flask

from mdvtools.auth.authutils import admin_required


def create_app(enable_auth: bool) -> Flask:
    app = Flask(__name__)
    app.secret_key = "test"
    app.config.update(
        ENABLE_AUTH=enable_auth,
        LOGIN_REDIRECT_URL="/login_dev",
    )

    @app.get("/protected-extension-route")
    @admin_required
    def protected_extension_route():
        return {"status": "allowed"}

    return app


def test_admin_required_allows_auth_disabled_development() -> None:
    response = (
        create_app(enable_auth=False).test_client().get("/protected-extension-route")
    )

    assert response.status_code == 200
    assert response.get_json() == {"status": "allowed"}


def test_admin_required_redirects_anonymous_user_to_mdv_login() -> None:
    response = (
        create_app(enable_auth=True).test_client().get("/protected-extension-route")
    )

    assert response.status_code == 302
    assert response.headers["Location"].endswith("/login_dev")


def test_admin_required_forbids_authenticated_non_admin() -> None:
    client = create_app(enable_auth=True).test_client()
    with client.session_transaction() as flask_session:
        flask_session["user"] = {"id": 1, "is_admin": False}

    assert client.get("/protected-extension-route").status_code == 403


def test_admin_required_allows_authenticated_admin() -> None:
    client = create_app(enable_auth=True).test_client()
    with client.session_transaction() as flask_session:
        flask_session["user"] = {"id": 2, "is_admin": True}

    response = client.get("/protected-extension-route")

    assert response.status_code == 200
    assert response.get_json() == {"status": "allowed"}
