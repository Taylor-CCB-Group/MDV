import os

import pandas as pd
import pytest

from mdvtools.mdvproject import MDVProject


@pytest.fixture
def project(tmp_path):
    """A throwaway project with one datasource of three columns: a, b, c.

    Uses pytest's tmp_path rather than tests/temp so each test is hermetic and
    nothing is written inside the repo.
    """
    project = MDVProject(os.path.join(str(tmp_path), "project"), delete_existing=True)
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
    project.add_datasource("test", df)
    return project


def column(project, field):
    return next(
        x for x in project.get_datasource_metadata("test")["columns"]
        if x.get("field") == field
    )


def test_renames_display_name_and_leaves_field_untouched(project):

    assert project.rename_column("test", "a", "Alpha") is True

    col = column(project, "a")
    assert col["name"] == "Alpha"
    assert col["field"] == "a"
    # other columns are untouched
    assert [x["name"] for x in project.get_datasource_metadata("test")["columns"]] == [
        "Alpha",
        "b",
        "c",
    ]


def test_rename_persists_to_disk(project):
    project.rename_column("test", "b", "Beta")

    # a fresh project object reads datasources.json again
    reloaded = MDVProject(project.dir)
    assert column(reloaded, "b")["name"] == "Beta"


def test_rename_strips_surrounding_whitespace(project):

    assert project.rename_column("test", "a", "  Alpha  ") is True
    assert column(project, "a")["name"] == "Alpha"


def test_rename_to_current_name_is_a_no_op(project):

    assert project.rename_column("test", "a", "Alpha") is True
    assert project.rename_column("test", "a", "Alpha") is False
    assert column(project, "a")["name"] == "Alpha"


def test_rename_preserves_other_column_metadata(project):
    before = dict(column(project, "a"))

    project.rename_column("test", "a", "Alpha")

    after = column(project, "a")
    assert set(after) == set(before)
    for key in before:
        if key != "name":
            assert after[key] == before[key], key


@pytest.mark.parametrize("new_name", ["", "   "])
def test_rename_rejects_empty_name(project, new_name):

    with pytest.raises(ValueError, match="required"):
        project.rename_column("test", "a", new_name)
    assert column(project, "a")["name"] == "a"


@pytest.mark.parametrize("new_name", ["b", "B", " b "])
def test_rename_rejects_duplicate_name_case_insensitively(project, new_name):

    with pytest.raises(ValueError, match="already exists"):
        project.rename_column("test", "a", new_name)
    assert column(project, "a")["name"] == "a"


def test_rename_allows_reusing_the_name_of_a_deleted_column(project):
    project.set_column_metadata("test", "c", "deleted", True)

    assert project.rename_column("test", "a", "c") is True
    assert column(project, "a")["name"] == "c"


def test_rename_rejects_unknown_column(project):

    with pytest.raises(AttributeError, match="not found"):
        project.rename_column("test", "nope", "Alpha")


def test_rename_rejects_unknown_datasource(project):

    with pytest.raises(AttributeError, match="datasource not found"):
        project.rename_column("nope", "a", "Alpha")


def test_rename_rejects_deleted_column(project):
    project.set_column_metadata("test", "a", "deleted", True)

    with pytest.raises(AttributeError, match="deleted"):
        project.rename_column("test", "a", "Alpha")


@pytest.mark.parametrize("spatial_key", ["sgindex", "sgtype"])
def test_rename_rejects_spatial_columns(project, spatial_key):
    project.set_column_metadata("test", "a", spatial_key, 0)

    with pytest.raises(AttributeError, match="spatial"):
        project.rename_column("test", "a", "Alpha")


def test_rename_rejects_legacy_column_with_no_field(project):
    ds = project.get_datasource_metadata("test")
    ds["columns"].append({"name": "legacy", "datatype": "text"})
    project.set_datasource_metadata(ds)

    with pytest.raises(AttributeError, match="no 'field' identifier"):
        project.rename_column("test", "legacy", "Legacy")


def test_rename_works_alongside_a_legacy_column_with_no_field(project):
    # a field-less column used to make the metadata lookups raise KeyError
    ds = project.get_datasource_metadata("test")
    ds["columns"].append({"name": "legacy", "datatype": "text"})
    project.set_datasource_metadata(ds)

    assert project.rename_column("test", "a", "Alpha") is True
    assert column(project, "a")["name"] == "Alpha"
