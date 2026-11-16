import copy
import hashlib
import os

import pandas as pd
import pytest
from flask import Flask

import mdvtools.mdvproject as mdvproject
from mdvtools.file_processing import ValidationError, validate_datasource
from mdvtools.mdvproject import MDVProject, datasource_display_name
from mdvtools.project_router import ProjectBlueprint
from mdvtools.server import create_app
from mdvtools.server_extension import MDVServerOptions
from mdvtools.websocket import mdv_socketio

UPLOAD = "upload.csv"


@pytest.fixture
def project(tmp_path):
    """A project with two datasources: `cells`, and `upload.csv` named after its file.

    `cells` links to `upload.csv` and the default view shows both, so there are
    references to the identifier that a rename must leave alone.
    """
    project = MDVProject(os.path.join(str(tmp_path), "project"), delete_existing=True)
    project.add_datasource("cells", pd.DataFrame({"x": [1, 2, 3]}))
    project.add_datasource(UPLOAD, pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}))

    cells = project.get_datasource_metadata("cells")
    cells["links"] = {UPLOAD: {"access_data": True, "index": "x"}}
    project.set_datasource_metadata(cells)
    return project


def entry(project, name):
    return project.get_datasource_metadata(name)


def h5_hash(project):
    with open(project.h5file, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


@pytest.fixture
def writes(monkeypatch):
    """Counts the JSON files the project writes."""
    written = []
    save_json = mdvproject.save_json

    def spy(file, data, *args, **kwargs):
        written.append(os.path.basename(file))
        return save_json(file, data, *args, **kwargs)

    monkeypatch.setattr(mdvproject, "save_json", spy)
    return written


# --- the rename itself ---


def test_sets_display_name_and_leaves_everything_else_untouched(project):
    datasources_before = copy.deepcopy(project.datasources)
    views_before = copy.deepcopy(project.views)
    h5_before = h5_hash(project)

    assert project.rename_datasource(UPLOAD, "Samples") is True

    datasources = project.datasources
    # only the renamed entry gains a key; the order, the other datasource and every
    # other key (including the identifier) are as they were
    assert [x["name"] for x in datasources] == ["cells", UPLOAD]
    assert datasources[0] == datasources_before[0]
    assert datasources[1] == {**datasources_before[1], "display_name": "Samples"}
    assert project.views == views_before
    assert h5_hash(project) == h5_before


def test_writes_only_datasources_json(project, writes):
    project.rename_datasource(UPLOAD, "Samples")

    assert writes == ["datasources.json"]


def test_persists_to_disk(project):
    project.rename_datasource(UPLOAD, "Samples")

    # a fresh project object reads datasources.json again
    assert entry(MDVProject(project.dir), UPLOAD)["display_name"] == "Samples"


def test_data_is_still_served_under_the_identifier(project):
    before = project.get_byte_data([{"field": "a"}], UPLOAD)

    project.rename_datasource(UPLOAD, "Samples")

    assert project.get_byte_data([{"field": "a"}], UPLOAD) == before
    with pytest.raises(AttributeError):
        project.get_datasource_metadata("Samples")


def test_strips_surrounding_whitespace(project):
    assert project.rename_datasource(UPLOAD, "  Samples\t") is True
    assert entry(project, UPLOAD)["display_name"] == "Samples"


def test_renaming_again_replaces_the_display_name(project):
    project.rename_datasource(UPLOAD, "Samples")

    assert project.rename_datasource(UPLOAD, "Patients") is True
    assert entry(project, UPLOAD)["display_name"] == "Patients"


def test_renaming_to_the_current_display_name_writes_nothing(project, writes):
    project.rename_datasource(UPLOAD, "Samples")
    writes.clear()

    assert project.rename_datasource(UPLOAD, " Samples ") is False
    assert writes == []


def test_renaming_to_its_own_name_when_never_renamed_writes_nothing(project, writes):
    assert project.rename_datasource(UPLOAD, UPLOAD) is False
    assert "display_name" not in entry(project, UPLOAD)
    assert writes == []


def test_renaming_back_to_its_own_name_keeps_the_field(project):
    project.rename_datasource(UPLOAD, "Samples")

    assert project.rename_datasource(UPLOAD, UPLOAD) is True
    assert entry(project, UPLOAD)["display_name"] == UPLOAD


@pytest.mark.parametrize("new_name", ["Upload.csv", "UPLOAD.CSV"])
def test_allows_changing_the_case_of_its_own_name(project, new_name):
    assert project.rename_datasource(UPLOAD, new_name) is True
    assert entry(project, UPLOAD)["display_name"] == new_name


def test_allows_changing_the_case_of_its_own_display_name(project):
    project.rename_datasource(UPLOAD, "Samples")

    assert project.rename_datasource(UPLOAD, "SAMPLES") is True
    assert entry(project, UPLOAD)["display_name"] == "SAMPLES"


def test_allows_a_display_name_another_datasource_gave_up(project):
    project.rename_datasource("cells", "Samples")
    project.rename_datasource("cells", "Cells")

    assert project.rename_datasource(UPLOAD, "Samples") is True


def test_keeps_unicode_names(project):
    assert project.rename_datasource(UPLOAD, "Échantillons 2026 – lot B") is True
    assert entry(MDVProject(project.dir), UPLOAD)["display_name"] == "Échantillons 2026 – lot B"


# --- rejected renames, which must leave the project as it was ---


@pytest.mark.parametrize("new_name", ["", "   ", "\t\n", None, 3, ["Samples"], {"name": "x"}])
def test_rejects_an_empty_or_non_text_name(project, writes, new_name):
    with pytest.raises(ValueError, match="Datasource name is required"):
        project.rename_datasource(UPLOAD, new_name)
    assert "display_name" not in entry(project, UPLOAD)
    assert writes == []


@pytest.mark.parametrize("datasource", ["missing", "Upload.csv", None])
def test_rejects_an_unknown_datasource(project, writes, datasource):
    # the identifier is matched exactly, as everywhere else in the project
    with pytest.raises(AttributeError, match="datasource not found"):
        project.rename_datasource(datasource, "Samples")
    assert writes == []


def test_rejects_a_datasource_given_by_its_display_name(project):
    project.rename_datasource(UPLOAD, "Samples")

    with pytest.raises(AttributeError, match="datasource not found"):
        project.rename_datasource("Samples", "Patients")


@pytest.mark.parametrize("new_name", ["cells", "CELLS", "  Cells  "])
def test_rejects_another_datasources_name_in_any_case(project, writes, new_name):
    with pytest.raises(ValueError, match="already exists"):
        project.rename_datasource(UPLOAD, new_name)
    assert "display_name" not in entry(project, UPLOAD)
    assert writes == []


@pytest.mark.parametrize("new_name", ["Samples", "samples", " SAMPLES "])
def test_rejects_another_datasources_display_name_in_any_case(project, new_name):
    project.rename_datasource("cells", "Samples")

    with pytest.raises(ValueError, match="already exists"):
        project.rename_datasource(UPLOAD, new_name)


def test_keeps_a_renamed_datasources_name_reserved(project):
    # otherwise `upload.csv` could be shown as "cells" while another datasource is cells
    project.rename_datasource("cells", "Cell table")

    with pytest.raises(ValueError, match="already exists"):
        project.rename_datasource(UPLOAD, "cells")


def test_ignores_a_corrupt_display_name_on_another_datasource(project):
    cells = entry(project, "cells")
    cells["display_name"] = 42
    project.set_datasource_metadata(cells)

    assert project.rename_datasource(UPLOAD, "42") is True


# --- the shown name, for entries written before display names existed ---


@pytest.mark.parametrize(
    "datasource, shown",
    [
        ({"name": "cells"}, "cells"),
        ({"name": "cells", "display_name": "Cells"}, "Cells"),
        ({"name": "cells", "display_name": "  Cells "}, "Cells"),
        ({"name": "cells", "display_name": ""}, "cells"),
        ({"name": "cells", "display_name": "   "}, "cells"),
        ({"name": "cells", "display_name": None}, "cells"),
        ({"name": "cells", "display_name": 42}, "cells"),
    ],
)
def test_datasource_display_name_falls_back_to_name(datasource, shown):
    assert datasource_display_name(datasource) == shown


def test_a_corrupt_display_name_on_the_datasource_itself_can_be_replaced(project):
    upload = entry(project, UPLOAD)
    upload["display_name"] = "   "
    project.set_datasource_metadata(upload)

    # it is shown as its name, so renaming to the name changes nothing
    assert project.rename_datasource(UPLOAD, UPLOAD) is False
    assert project.rename_datasource(UPLOAD, "Samples") is True
    assert entry(project, UPLOAD)["display_name"] == "Samples"


# --- the route ---


@pytest.fixture
def client(project):
    app = Flask(__name__)
    mdv_socketio(app)
    ProjectBlueprint.register_app(app)
    create_app(project, MDVServerOptions(app=app, open_browser=False, websocket=False))
    yield app.test_client()
    ProjectBlueprint.blueprints.pop(project.id, None)


def url(project):
    return f"/project/{project.id}/rename_datasource"


def post(client, url, **kwargs):
    # the project route is registered with a trailing slash, so Flask answers 308 first
    return client.post(url, follow_redirects=True, **kwargs)


def test_route_renames(client, project):
    response = post(client, url(project), json={"datasource": UPLOAD, "name": " Samples "})

    assert response.status_code == 200
    assert response.json == {"success": True, "changed": True, "display_name": "Samples"}
    assert entry(project, UPLOAD)["display_name"] == "Samples"


def test_route_reports_an_unchanged_name(client, project):
    post(client, url(project), json={"datasource": UPLOAD, "name": "Samples"})

    response = post(client, url(project), json={"datasource": UPLOAD, "name": "Samples"})

    assert response.json == {"success": True, "changed": False, "display_name": "Samples"}


def test_route_still_serves_data_under_the_identifier(client, project):
    before = project.get_byte_data([{"field": "a"}], UPLOAD)
    post(client, url(project), json={"datasource": UPLOAD, "name": "Samples"})

    data_url = f"/project/{project.id}/get_data"
    data = post(client, data_url, json={"columns": [{"field": "a"}], "data_source": UPLOAD})

    assert data.status_code == 200
    assert data.data == before


@pytest.mark.parametrize(
    "kwargs",
    [
        {"json": {"datasource": UPLOAD}},
        {"json": {"name": "Samples"}},
        {"json": {}},
        {"json": ["datasource", "name"]},
        {"data": "datasource=upload.csv&name=Samples", "content_type": "text/plain"},
        {"data": "{not json", "content_type": "application/json"},
    ],
)
def test_route_rejects_a_malformed_body(client, project, kwargs):
    response = post(client, url(project), **kwargs)

    assert response.status_code == 400
    assert response.json == {
        "success": False,
        "error": "Request must contain 'datasource' and 'name'",
    }


@pytest.mark.parametrize(
    "body, error",
    [
        ({"datasource": UPLOAD, "name": "  "}, "Datasource name is required"),
        ({"datasource": UPLOAD, "name": None}, "Datasource name is required"),
        ({"datasource": UPLOAD, "name": "Cells"}, "A datasource called 'Cells' already exists"),
        ({"datasource": "missing", "name": "Samples"}, "missing datasource not found"),
    ],
)
def test_route_reports_errors(client, project, body, error):
    response = post(client, url(project), json=body)

    assert response.status_code == 400
    assert response.json == {"success": False, "error": error}
    assert "display_name" not in entry(project, UPLOAD)


def test_route_requires_an_editable_project(client, project):
    # the auth hooks that enforce access levels are not installed in tests, so check that
    # the route asks for one
    routes = ProjectBlueprint.blueprints[project.id].routes
    options = [opts for rule, (_, opts) in routes.items() if rule.match("/rename_datasource")]

    assert options == [{"access_level": "editable", "methods": ["POST"]}]


# --- uploads, which must not add a second datasource shown under the same name ---


def upload_request(name, replace=False):
    return {"name": name, "replace": replace}


@pytest.fixture
def editable(project):
    project.set_editable(True)
    return project


@pytest.mark.parametrize("name", ["Samples", "samples", " SAMPLES "])
def test_upload_rejects_a_name_a_datasource_is_shown_as(editable, name):
    editable.rename_datasource(UPLOAD, "Samples")

    with pytest.raises(ValidationError, match="already shown as"):
        validate_datasource(editable, upload_request(name))
    with pytest.raises(ValidationError, match="already shown as"):
        validate_datasource(editable, upload_request(name, replace=True))


def test_upload_can_still_replace_a_renamed_datasource_by_identifier(editable):
    editable.rename_datasource(UPLOAD, "Samples")

    assert validate_datasource(editable, upload_request(UPLOAD, replace=True))["name"] == UPLOAD
    with pytest.raises(ValidationError, match="already exists"):
        validate_datasource(editable, upload_request(UPLOAD))


def test_upload_rules_are_unchanged_for_datasources_never_renamed(editable):
    # a name differing only in case from an identifier was allowed before display names
    assert validate_datasource(editable, upload_request("CELLS"))["name"] == "CELLS"
    assert validate_datasource(editable, upload_request("new.csv"))["name"] == "new.csv"
