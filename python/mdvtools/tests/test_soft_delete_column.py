import os

import pandas as pd
import pytest

from mdvtools.mdvproject import MDVProject


@pytest.fixture
def project(tmp_path):
    """A throwaway project with one datasource of three columns: a, b, c, and no views.

    `add_to_view=None` matters: the default (`"default"`) makes add_datasource build a
    view containing a table plot whose params are every column's field, which would
    block every deletion. A curation script has to do the same - curate first, then
    generate views. See test_a_generated_default_view_blocks_deletion below.

    Uses pytest's tmp_path rather than tests/temp so each test is hermetic and
    nothing is written inside the repo.
    """
    project = MDVProject(os.path.join(str(tmp_path), "project"), delete_existing=True)
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
    project.add_datasource("test", df, add_to_view=None)
    return project


def column(project, field):
    return next(
        x for x in project.get_datasource_metadata("test")["columns"]
        if x.get("field") == field
    )


def add_view(project, name, charts, datasource="test"):
    project.set_view(name, {"initialCharts": {datasource: charts}})


# --- tombstoning -------------------------------------------------------------


def test_soft_delete_sets_the_tombstone(project):
    assert project.soft_delete_column("test", "a") is True

    col = column(project, "a")
    assert col["deleted"] is True
    # the column and its identity are still there - this is a tombstone, not a removal
    assert col["field"] == "a"
    assert col["name"] == "a"


def test_soft_delete_keeps_the_column_in_the_datasource(project):
    before = len(project.get_datasource_metadata("test")["columns"])

    project.soft_delete_column("test", "a")

    after = project.get_datasource_metadata("test")["columns"]
    assert len(after) == before
    assert [x["field"] for x in after] == ["a", "b", "c"]


def test_soft_delete_leaves_other_columns_untouched(project):
    project.soft_delete_column("test", "a")

    for field in ("b", "c"):
        assert "deleted" not in column(project, field)


def test_soft_delete_preserves_other_column_metadata(project):
    before = dict(column(project, "a"))

    project.soft_delete_column("test", "a")

    after = column(project, "a")
    for key in before:
        assert after[key] == before[key], key


def test_soft_delete_persists_to_disk(project):
    project.soft_delete_column("test", "a")

    reloaded = MDVProject(project.dir)
    assert column(reloaded, "a")["deleted"] is True


def test_soft_delete_of_an_already_deleted_column_does_nothing(project):
    assert project.soft_delete_column("test", "a") is True
    assert project.soft_delete_column("test", "a") is False
    assert column(project, "a")["deleted"] is True


# --- rejections --------------------------------------------------------------


def test_soft_delete_rejects_unknown_column(project):
    with pytest.raises(AttributeError, match="not found"):
        project.soft_delete_column("test", "nope")


def test_soft_delete_rejects_unknown_datasource(project):
    with pytest.raises(AttributeError, match="datasource not found"):
        project.soft_delete_column("nope", "a")


@pytest.mark.parametrize("spatial_key", ["sgindex", "sgtype"])
def test_soft_delete_rejects_spatial_columns(project, spatial_key):
    project.set_column_metadata("test", "a", spatial_key, 0)

    with pytest.raises(AttributeError, match="spatial"):
        project.soft_delete_column("test", "a")
    assert "deleted" not in column(project, "a")


def test_soft_delete_rejects_legacy_column_with_no_field(project):
    ds = project.get_datasource_metadata("test")
    ds["columns"].append({"name": "legacy", "datatype": "text"})
    project.set_datasource_metadata(ds)

    with pytest.raises(AttributeError, match="no 'field' identifier"):
        project.soft_delete_column("test", "legacy")


# --- blocked by views --------------------------------------------------------


def test_soft_delete_is_blocked_by_a_chart_param(project):
    add_view(project, "Overview", [{"type": "table_chart", "param": ["a", "b"]}])

    with pytest.raises(ValueError, match="blocked"):
        project.soft_delete_column("test", "a")
    assert "deleted" not in column(project, "a")


def test_blocked_message_names_the_view_and_the_chart(project):
    add_view(project, "Overview", [
        {"type": "table_chart", "title": "Cell table", "param": ["a"]},
    ])

    with pytest.raises(ValueError) as excinfo:
        project.soft_delete_column("test", "a")

    message = str(excinfo.value)
    assert "still used by other charts or views" in message
    assert "Overview" in message
    assert "Cell table" in message


def test_soft_delete_is_blocked_by_a_nested_setting(project):
    # proves the scan is recursive rather than looking only at `param`
    add_view(project, "Scatter", [
        {"type": "wgl_scatter_plot", "param": ["b", "c"], "color_by": {"column": "a"}},
    ])

    with pytest.raises(ValueError, match="blocked"):
        project.soft_delete_column("test", "a")


def test_soft_delete_is_blocked_by_any_of_several_views(project):
    add_view(project, "First", [{"type": "table_chart", "param": ["b"]}])
    add_view(project, "Second", [{"type": "table_chart", "param": ["a"]}])

    with pytest.raises(ValueError, match="Second"):
        project.soft_delete_column("test", "a")


def test_soft_delete_succeeds_when_views_use_other_columns(project):
    add_view(project, "Overview", [{"type": "table_chart", "param": ["b", "c"]}])

    assert project.soft_delete_column("test", "a") is True


def test_soft_delete_ignores_charts_on_another_datasource(project):
    add_view(project, "Other", [{"type": "table_chart", "param": ["a"]}], datasource="other")

    assert project.soft_delete_column("test", "a") is True


def test_soft_delete_succeeds_when_there_are_no_views(project):
    assert project.soft_delete_column("test", "a") is True


def test_scan_matches_whole_strings_only(project):
    # "a" must not match "abc" - the scan compares strings exactly
    add_view(project, "Overview", [{"type": "table_chart", "param": ["abc"]}])

    assert project.soft_delete_column("test", "a") is True


def test_scan_over_reports_rather_than_missing_a_reference(project):
    # Documented fail-closed behaviour: a field name appearing in an unrelated string
    # (here a chart title) blocks the delete. Over-blocking is recoverable; a missed
    # reference is not. See ADR-0003.
    add_view(project, "Overview", [{"type": "table_chart", "title": "a", "param": ["b"]}])

    with pytest.raises(ValueError, match="blocked"):
        project.soft_delete_column("test", "a")


def test_a_generated_default_view_blocks_deletion(tmp_path):
    """add_datasource's default view references every column, so it blocks everything.

    This is the real-world trap: a script that creates a datasource with the default
    add_to_view="default" cannot then hide any of its columns.
    """
    project = MDVProject(os.path.join(str(tmp_path), "blocked"), delete_existing=True)
    project.add_datasource("test", pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}))

    with pytest.raises(ValueError, match="blocked"):
        project.soft_delete_column("test", "a")


def test_generated_views_exclude_tombstoned_columns(project):
    """create_view_with_all_datasources must not list a hidden column in its table plot.

    Otherwise the curate-then-generate order still produces a view referencing columns
    that MDV will not render, and the next deletion is blocked by that generated view.
    """
    project.soft_delete_column("test", "a")

    view_name = project.create_view_with_all_datasources("Generated")

    charts = project.get_view(view_name)["initialCharts"]["test"]
    assert len(charts) == 1
    # key-agnostic: nothing in the generated view references the hidden column,
    # while the surviving columns are referenced
    assert project._find_column_references("test", "a") == []
    assert project._find_column_references("test", "b") != []
    assert project._find_column_references("test", "c") != []
