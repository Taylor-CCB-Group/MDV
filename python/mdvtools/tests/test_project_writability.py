import logging
import os

from mdvtools.mdvproject import MDVProject
from mdvtools.server import _apply_project_writability


def test_unwritable_paths_names_every_failed_access_check(tmp_path, monkeypatch):
    project = MDVProject(str(tmp_path / "project"))
    with open(project.h5file, "wb"):
        pass

    blocked_paths = {
        project.statefile,
        project.dir,
        project.viewsfile,
        project.h5file,
    }
    monkeypatch.setattr(
        os,
        "access",
        lambda path, _mode: path not in blocked_paths,
    )

    assert project.unwritable_paths == [
        project.statefile,
        project.dir,
        project.viewsfile,
        project.h5file,
    ]
    assert not project.writable


def test_set_editable_warns_with_failed_paths(tmp_path, monkeypatch, caplog):
    project = MDVProject(str(tmp_path / "project"))
    monkeypatch.setattr(
        MDVProject,
        "unwritable_paths",
        property(lambda self: [self.statefile, self.viewsfile]),
    )
    caplog.set_level(logging.WARNING, logger="mdvtools.mdvproject")

    project.set_editable()

    assert "Cannot set project" in caplog.text
    assert project.statefile in caplog.text
    assert project.viewsfile in caplog.text


def test_state_writability_warns_once_and_serves_view_permission(
    tmp_path,
    monkeypatch,
    caplog,
):
    project = MDVProject(str(tmp_path / "project"), id="7")
    monkeypatch.setattr(
        MDVProject,
        "unwritable_paths",
        property(lambda self: [self.statefile]),
    )
    caplog.set_level(logging.WARNING, logger="mdvtools.server")
    first_state = {"permission": "edit"}
    second_state = {"permission": "edit"}

    _apply_project_writability(project, first_state, websocket_enabled=True)
    _apply_project_writability(project, second_state, websocket_enabled=True)

    assert first_state == {"permission": "view", "websocket": False}
    assert second_state == {"permission": "view", "websocket": False}
    warnings = [
        record
        for record in caplog.records
        if "Serving project read-only" in record.getMessage()
    ]
    assert len(warnings) == 1
    assert project.statefile in warnings[0].getMessage()
