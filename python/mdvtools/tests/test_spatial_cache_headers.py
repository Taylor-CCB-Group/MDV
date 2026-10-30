"""Tests for the caching headers on files served from a project's spatial folder."""

import os

import pytest
from flask import Flask

from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprint
from mdvtools.server import create_app
from mdvtools.server_extension import MDVServerOptions
from mdvtools.websocket import mdv_socketio


@pytest.fixture
def client_and_project(tmp_path):
    project = MDVProject(str(tmp_path / "project"))
    chunk_dir = os.path.join(project.dir, "spatial", "sample.zarr", "images", "img", "0", "c", "0", "0")
    os.makedirs(chunk_dir)
    with open(os.path.join(chunk_dir, "0"), "wb") as f:
        f.write(b"chunk bytes")

    app = Flask(__name__)
    mdv_socketio(app)
    ProjectBlueprint.register_app(app)
    create_app(project, MDVServerOptions(app=app, open_browser=False, websocket=False))
    yield app.test_client(), project
    ProjectBlueprint.blueprints.pop(project.id, None)


def _get_chunk(client, project, etag=None):
    headers = {"If-None-Match": etag} if etag else {}
    return client.get(
        f"/project/{project.id}/spatial/sample.zarr/images/img/0/c/0/0/0",
        headers=headers,
        follow_redirects=True,
    )


def test_spatial_files_are_cached_privately_for_an_hour(client_and_project):
    client, project = client_and_project

    response = _get_chunk(client, project)

    assert response.status_code == 200
    assert response.data == b"chunk bytes"
    assert response.cache_control.max_age == 3600
    assert response.cache_control.private is True
    assert response.cache_control.public is False
    assert response.cache_control.no_cache is None
    assert _get_chunk(client, project, response.headers["ETag"]).status_code == 304


def test_spatial_files_revalidate_every_request_when_auth_is_on(client_and_project):
    client, project = client_and_project
    client.application.config["ENABLE_AUTH"] = True

    response = _get_chunk(client, project)

    assert response.status_code == 200
    assert response.cache_control.max_age is None
    assert response.cache_control.private is True
    assert response.cache_control.public is False
    assert response.cache_control.no_cache is True
    assert _get_chunk(client, project, response.headers["ETag"]).status_code == 304


def test_missing_spatial_files_are_not_cached(client_and_project):
    client, project = client_and_project

    response = client.get(f"/project/{project.id}/spatial/sample.zarr/missing", follow_redirects=True)

    assert response.status_code == 404
    assert response.cache_control.max_age is None
