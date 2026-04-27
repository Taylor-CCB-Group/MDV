"""Unit tests for project_pre_dispatch_checks (project_auth module)."""

import pytest
from flask import Flask
from unittest.mock import patch, MagicMock

from mdvtools.auth.project_auth import project_pre_dispatch_checks


@pytest.fixture
def app():
    """Create and configure a new app instance for each test."""
    app = Flask("test_project_auth")
    app.config["TESTING"] = True
    app.config["ENABLE_AUTH"] = True
    app.secret_key = "test_secret"
    return app


@pytest.fixture
def mock_project():
    """Project model instance with access_level and accessed_timestamp."""
    p = MagicMock()
    p.access_level = "editable"
    p.accessed_timestamp = None
    return p


def test_project_not_found_returns_404(app, mock_project):
    """Project not found returns 404 before any user lookup."""
    with patch("mdvtools.auth.project_auth.ProjectService") as mock_ps:
        mock_ps.get_project_by_id.return_value = None
        with app.app_context():
            with app.test_request_context():
                rv = project_pre_dispatch_checks(project_id="999", options={})
        assert rv is not None
        json_data, status_code = rv
        assert status_code == 404
        mock_ps.get_project_by_id.assert_called_once_with("999")


def test_auth_disabled_allows_without_user(app, mock_project):
    """When ENABLE_AUTH is False, request is allowed without user in session."""
    with patch("mdvtools.auth.project_auth.ProjectService") as mock_ps:
        mock_ps.get_project_by_id.return_value = mock_project
        mock_ps.set_project_accessed_timestamp.return_value = None
        with app.app_context():
            app.config["ENABLE_AUTH"] = False
            with app.test_request_context():
                # No session['user']
                rv = project_pre_dispatch_checks(project_id="1", options={})
        assert rv is None


def test_auth_enabled_no_user_returns_401(app, mock_project):
    """When ENABLE_AUTH is True and no user in session, returns 401."""
    with patch("mdvtools.auth.project_auth.ProjectService") as mock_ps:
        mock_ps.get_project_by_id.return_value = mock_project
        mock_ps.set_project_accessed_timestamp.return_value = None
        with app.app_context():
            app.config["ENABLE_AUTH"] = True
            with app.test_request_context():
                # session has no 'user'
                rv = project_pre_dispatch_checks(project_id="1", options={})
        assert rv is not None
        json_data, status_code = rv
        assert status_code == 401


def test_auth_enabled_user_with_no_project_permission_returns_403(app, mock_project):
    """When ENABLE_AUTH is True and user has no permission for project, returns 403."""
    with patch("mdvtools.auth.project_auth.ProjectService") as mock_ps, \
         patch("mdvtools.auth.project_auth.user_project_cache", {42: {}}) as cache:
        mock_ps.get_project_by_id.return_value = mock_project
        mock_ps.set_project_accessed_timestamp.return_value = None
        with app.app_context():
            app.config["ENABLE_AUTH"] = True
            with app.test_request_context():
                from flask import session
                session["user"] = {"id": 42}
                # user 42 has no entry for project 1 in cache (only empty dict for user 42)
                rv = project_pre_dispatch_checks(project_id="1", options={})
        assert rv is not None
        json_data, status_code = rv
        assert status_code == 403
        assert b"You do not have access to this project" in json_data.get_data()


def test_auth_enabled_user_with_no_read_flags_returns_403(app, mock_project):
    """User has project entry but can_read/can_write/is_owner all False -> 403."""
    with patch("mdvtools.auth.project_auth.ProjectService") as mock_ps, \
         patch("mdvtools.auth.project_auth.user_project_cache") as cache:
        cache.get.return_value = {
            1: {"can_read": False, "can_write": False, "is_owner": False},
        }
        mock_ps.get_project_by_id.return_value = mock_project
        mock_ps.set_project_accessed_timestamp.return_value = None
        with app.app_context():
            app.config["ENABLE_AUTH"] = True
            with app.test_request_context():
                from flask import session
                session["user"] = {"id": 99}
                cache.get.return_value = {1: {"can_read": False, "can_write": False, "is_owner": False}}
                rv = project_pre_dispatch_checks(project_id="1", options={})
        assert rv is not None
        json_data, status_code = rv
        assert status_code == 403


def test_auth_enabled_can_read_only_allowed_for_non_editable_route(app, mock_project):
    """User with can_read only is allowed for route without access_level='editable'."""
    with patch("mdvtools.auth.project_auth.ProjectService") as mock_ps, \
         patch("mdvtools.auth.project_auth.user_project_cache") as cache:
        cache.get.return_value = {
            1: {"can_read": True, "can_write": False, "is_owner": False},
        }
        mock_ps.get_project_by_id.return_value = mock_project
        mock_ps.set_project_accessed_timestamp.return_value = None
        with app.app_context():
            app.config["ENABLE_AUTH"] = True
            with app.test_request_context():
                from flask import session
                session["user"] = {"id": 10}
                rv = project_pre_dispatch_checks(project_id="1", options={})
        assert rv is None


def test_auth_enabled_can_read_only_denied_for_editable_route(app, mock_project):
    """User with can_read only gets 403 for access_level='editable' route."""
    with patch("mdvtools.auth.project_auth.ProjectService") as mock_ps, \
         patch("mdvtools.auth.project_auth.user_project_cache") as cache:
        cache.get.return_value = {
            1: {"can_read": True, "can_write": False, "is_owner": False},
        }
        mock_ps.get_project_by_id.return_value = mock_project
        mock_ps.set_project_accessed_timestamp.return_value = None
        with app.app_context():
            app.config["ENABLE_AUTH"] = True
            with app.test_request_context():
                from flask import session
                session["user"] = {"id": 10}
                rv = project_pre_dispatch_checks(project_id="1", options={"access_level": "editable"})
        assert rv is not None
        json_data, status_code = rv
        assert status_code == 403


def test_auth_enabled_can_write_allowed_for_editable_route(app, mock_project):
    """User with can_write is allowed for access_level='editable' route."""
    with patch("mdvtools.auth.project_auth.ProjectService") as mock_ps, \
         patch("mdvtools.auth.project_auth.user_project_cache") as cache:
        cache.get.return_value = {
            1: {"can_read": True, "can_write": True, "is_owner": False},
        }
        mock_ps.get_project_by_id.return_value = mock_project
        mock_ps.set_project_accessed_timestamp.return_value = None
        mock_ps.set_project_update_timestamp.return_value = None
        with app.app_context():
            app.config["ENABLE_AUTH"] = True
            with app.test_request_context():
                from flask import session
                session["user"] = {"id": 10}
                rv = project_pre_dispatch_checks(project_id="1", options={"access_level": "editable"})
        assert rv is None
        mock_ps.set_project_update_timestamp.assert_called_once_with("1")


def test_auth_enabled_is_owner_allowed_for_editable_route(app, mock_project):
    """User with is_owner is allowed for access_level='editable' route."""
    with patch("mdvtools.auth.project_auth.ProjectService") as mock_ps, \
         patch("mdvtools.auth.project_auth.user_project_cache") as cache:
        cache.get.return_value = {
            1: {"can_read": True, "can_write": False, "is_owner": True},
        }
        mock_ps.get_project_by_id.return_value = mock_project
        mock_ps.set_project_accessed_timestamp.return_value = None
        mock_ps.set_project_update_timestamp.return_value = None
        with app.app_context():
            app.config["ENABLE_AUTH"] = True
            with app.test_request_context():
                from flask import session
                session["user"] = {"id": 10}
                rv = project_pre_dispatch_checks(project_id="1", options={"access_level": "editable"})
        assert rv is None
