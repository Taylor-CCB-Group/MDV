"""Tests for the caching headers on files served from a project's spatial folder."""

import os
from typing import BinaryIO

import pytest
from flask import Flask
from flask.testing import FlaskClient

from mdvtools import server_utils
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


class _NonSeekableFileWrapper:
    """Stands in for gunicorn's wsgi.file_wrapper, which has no seekable()."""

    bytes_read = 0

    def __init__(self, filelike: BinaryIO, blksize: int = 8192) -> None:
        self.filelike = filelike
        self.blksize = blksize
        self.close = filelike.close

    def __iter__(self) -> "_NonSeekableFileWrapper":
        return self

    def __next__(self) -> bytes:
        data = self.filelike.read(self.blksize)
        if not data:
            raise StopIteration
        type(self).bytes_read += len(data)
        return data


def test_spatial_range_requests_seek_instead_of_scanning_from_the_start(
    client_and_project: tuple[FlaskClient, MDVProject], monkeypatch: pytest.MonkeyPatch
) -> None:
    client, project = client_and_project
    content = bytes(range(256)) * 4096
    with open(os.path.join(project.dir, "spatial", "store.zarr.zip"), "wb") as f:
        f.write(content)
    _NonSeekableFileWrapper.bytes_read = 0
    served_bytes_read = 0

    class _CountingBlocks(server_utils._SeekableBlocks):
        def __next__(self) -> bytes:
            nonlocal served_bytes_read
            block = super().__next__()
            served_bytes_read += len(block)
            return block

    monkeypatch.setattr(server_utils, "_SeekableBlocks", _CountingBlocks)

    response = client.get(
        # trailing slash avoids the redirect, which would drop environ_overrides
        f"/project/{project.id}/spatial/store.zarr.zip/",
        headers={"Range": "bytes=900000-900999"},
        environ_overrides={"wsgi.file_wrapper": _NonSeekableFileWrapper},
    )

    assert response.status_code == 206
    assert response.data == content[900000:901000]
    assert response.headers["Content-Range"] == f"bytes 900000-900999/{len(content)}"
    assert _NonSeekableFileWrapper.bytes_read == 0
    assert 0 < served_bytes_read < 900000


def test_spatial_suffix_range_returns_the_end_of_the_file(client_and_project: tuple[FlaskClient, MDVProject]) -> None:
    client, project = client_and_project
    content = bytes(range(256)) * 100
    with open(os.path.join(project.dir, "spatial", "store.zarr.zip"), "wb") as f:
        f.write(content)

    response = client.get(
        f"/project/{project.id}/spatial/store.zarr.zip",
        headers={"Range": "bytes=-100"},
        follow_redirects=True,
    )

    assert response.status_code == 206
    assert response.data == content[-100:]
    assert response.cache_control.max_age == 3600


def test_spatial_unsatisfiable_range_is_rejected(client_and_project: tuple[FlaskClient, MDVProject]) -> None:
    client, project = client_and_project
    with open(os.path.join(project.dir, "spatial", "store.zarr.zip"), "wb") as f:
        f.write(b"small")

    response = client.get(
        f"/project/{project.id}/spatial/store.zarr.zip",
        headers={"Range": "bytes=1000-2000"},
        follow_redirects=True,
    )

    assert response.status_code == 416


def test_spatial_range_and_full_responses_share_an_etag(client_and_project: tuple[FlaskClient, MDVProject]) -> None:
    client, project = client_and_project
    with open(os.path.join(project.dir, "spatial", "store.zarr.zip"), "wb") as f:
        f.write(bytes(range(256)) * 100)
    url = f"/project/{project.id}/spatial/store.zarr.zip"

    full = client.get(url, follow_redirects=True)
    ranged = client.get(url, headers={"Range": "bytes=10-19"}, follow_redirects=True)

    assert ranged.status_code == 206
    assert ranged.headers["ETag"] == full.headers["ETag"]
    assert ranged.headers["Last-Modified"] == full.headers["Last-Modified"]
