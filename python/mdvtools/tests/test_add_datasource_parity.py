"""Parity tests for add_datasource (pandas) vs add_datasource_polars. 

Both methods accept either in-memory dataframe or a path to a delimited text file, and both
build a default TablePlot when add_to_view is set. These tests pin shared behaviour so a fix
landing on only one of the two implemntations is caught.
"""
import copy
import pandas as pd
import pytest

from mdvtools.mdvproject import MDVProject

def _write_tsv(tmp_path, filename: str = "data.tsv") -> str:
    """Helper to write a simple dataframe to a tsv file."""
    path = tmp_path / filename
    pd.DataFrame(
        {
            "keep_me": [1, 2, 3],
            "and_me": ["x", "y", "z"],
            "ignore_me": [10.0, 20.0, 30.0],
        }
    ).to_csv(path, sep="\t", index=False)
    return str(path)

def _table_plot_params(
        project: MDVProject, datasource_name: str, view_name: str = "Default View"
) -> list[str]:
    """Helper to get the datasource parameters of the default TablePlot created by add_datasource."""
    view = project.get_view(view_name)
    charts = view["initialCharts"][datasource_name]

    assert len(charts) == 1, (
        f"expected one chart in view {view_name} got {len(charts)}"
    )
    return list(charts[0]["param"])

SUPPLIED_COLUMNS = [
    {"name": "keep_me", "field": "keep_me", "datatype": "integer"},
    {"name": "and_me", "field": "and_me", "datatype": "text"},
]
EXPECTED_PARAMS = ["keep_me", "and_me"]

@pytest.mark.filterwarnings("ignore:Failed to add column")
def test_pandas_table_plot_params_respects_supplied_columns_only(tmp_path):
    """When supplied_columns_only=True, the default TablePlot's param list contains only 
    the supplied subset, not every column in the source file."""
    csv = _write_tsv(tmp_path)
    p = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    p.add_datasource(
        name="ds",
        dataframe=csv,
        columns=copy.deepcopy(SUPPLIED_COLUMNS),
        supplied_columns_only=True,
        add_to_view="Default View",
    )
    assert _table_plot_params(p, "ds", "Default View") == EXPECTED_PARAMS

@pytest.mark.filterwarnings("ignore:Failed to add column")
def test_polars_table_plot_params_respects_supplied_columns_only(tmp_path):
    """When supplied_columns_only=True, the default TablePlot's param list contains only 
    the supplied subset, not every column in the source file."""
    csv = _write_tsv(tmp_path)
    p = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    p.add_datasource_polars(
        name="ds",
        dataframe=csv,
        columns=copy.deepcopy(SUPPLIED_COLUMNS),
        supplied_columns_only=True,
        add_to_view="Default View",
    )
    assert _table_plot_params(p, "ds", "Default View") == EXPECTED_PARAMS