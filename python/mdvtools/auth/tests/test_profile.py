"""/profile answers from this deployment's users table, without calling the identity provider."""

import pytest
from flask import Flask

from mdvtools.auth.authutils import register_before_request_auth
from mdvtools.auth.register_auth_routes import register_auth_routes
from mdvtools.dbutils.dbmodels import User, db

USER = {"id": 7, "auth_id": "auth0|someone", "email": "someone@example.com", "is_admin": True}


@pytest.fixture
def app():
    app = Flask("test_profile")
    app.config.update(
        TESTING=True,
        ENABLE_AUTH=True,
        DEFAULT_AUTH_METHOD="auth0",
        LOGIN_REDIRECT_URL="/login_dev",
        SQLALCHEMY_DATABASE_URI="sqlite:///:memory:",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
    )
    app.secret_key = "test_secret"
    db.init_app(app)
    register_before_request_auth(app)
    register_auth_routes(app)

    with app.app_context():
        db.create_all()
        db.session.add(
            User(
                **USER,
                first_name="Ada",
                last_name="Lovelace",
                institution="University of Oxford",
                is_active=True,
                password="",
            )
        )
        db.session.commit()
        yield app
        db.session.remove()
        db.drop_all()


def test_profile_comes_from_the_database(app):
    client = app.test_client()
    with client.session_transaction() as sess:
        # An Auth0 access token that has expired. /profile must not need it.
        sess["user"] = dict(USER)
        sess["token"] = {"access_token": "expired"}

    rv = client.get("/profile")

    assert rv.status_code == 200
    assert rv.get_json() == {
        "name": "Ada Lovelace",
        "email": "someone@example.com",
        "association": "University of Oxford",
    }
