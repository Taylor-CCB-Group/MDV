"""Columns with pandas' nullable `string` dtype, as anndata reads a `nullable-string-array`
(for example the obs index of stores written by spatialdata-io's Atera reader)."""

import pandas as pd

from mdvtools.mdvproject import MDVProject


def test_string_dtype_columns_become_text_and_unique(tmp_path):
    mdv = MDVProject(str(tmp_path / "project"), delete_existing=True)
    df = pd.DataFrame(
        {
            "barcode": pd.array(["aaa-1", "aab-1", "aac-1", "aad-1"], dtype="string"),
            "kind": pd.array(["a", "b", None, "a"], dtype="string"),
        }
    )
    mdv.add_datasource("cells", df, [{"name": "barcode", "field": "barcode", "datatype": "unique"}])

    assert mdv.get_column_metadata("cells", "barcode")["datatype"] == "unique"
    assert mdv.get_column_metadata("cells", "kind")["datatype"] == "text"
    assert list(mdv.get_column("cells", "barcode")) == ["aaa-1", "aab-1", "aac-1", "aad-1"]
    assert list(mdv.get_column("cells", "kind")) == ["a", "b", "ND", "a"]
