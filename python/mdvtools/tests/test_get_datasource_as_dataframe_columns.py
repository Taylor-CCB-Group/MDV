"""Optional ``columns`` filter on get_datasource_as_dataframe (LLM-generated scripts)."""

import numpy
import pandas

from mdvtools.mdvproject import MDVProject


def test_get_datasource_as_dataframe_columns_subset(tmp_path):
    p = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    p.add_datasource("ds", pandas.DataFrame({"a": [1, 2], "b": [3, 4]}))
    full = p.get_datasource_as_dataframe("ds")
    assert list(full.columns) == ["a", "b"]
    sub = p.get_datasource_as_dataframe("ds", columns=["b"])
    assert list(sub.columns) == ["b"]
    assert sub["b"].tolist() == [3, 4]


def test_get_datasource_as_dataframe_columns_order(tmp_path):
    p = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    p.add_datasource("ds", pandas.DataFrame({"x": [1], "y": [2], "z": [3]}))
    sub = p.get_datasource_as_dataframe("ds", columns=["z", "x"])
    assert list(sub.columns) == ["z", "x"]


def test_get_datasource_as_dataframe_unknown_column_raises(tmp_path):
    p = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    p.add_datasource("ds", pandas.DataFrame({"a": [1]}))
    try:
        p.get_datasource_as_dataframe("ds", columns=["nope"])
    except AttributeError as e:
        assert "nope" in str(e)
    else:
        raise AssertionError("expected AttributeError")


def test_get_datasource_as_dataframe_chart_wrapper_plus_metadata(tmp_path):
    """Chart FieldName wrappers address rows-as-columns matrix columns, not metadata field ids."""
    p = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    n_cells, n_genes = 4, 3
    obs = pandas.DataFrame(
        {"cell_id": [f"c{i}" for i in range(n_cells)], "leiden": [0, 0, 1, 1]}
    )
    p.add_datasource("cells", obs, columns=[{"name": "cell_id", "datatype": "unique"}])
    var = pandas.DataFrame({"name": ["A", "B", "C"]})
    p.add_datasource("genes", var)
    p.add_rows_as_columns_link("cells", "genes", "name", "GE")
    x = numpy.arange(n_cells * n_genes, dtype=numpy.float32).reshape(n_cells, n_genes)
    p.add_rows_as_columns_subgroup("cells", "genes", "gs", x, name="scores", sparse=False)
    w = "gs|B(gs)|1"
    df = p.get_datasource_as_dataframe("cells", columns=["leiden", w])
    assert len(df) == n_cells
    numpy.testing.assert_array_almost_equal(df[w].to_numpy(), x[:, 1])
    assert df["leiden"].tolist() == [0, 0, 1, 1]
