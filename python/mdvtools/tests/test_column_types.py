"""The column-type contract (ADR-0024) as the Python package reads it."""

import numpy as np
import pandas as pd
import polars as pl
import pytest

from mdvtools.column_types import LEGACY_DATATYPES, NUMPY_STORAGE_DTYPES, infer_datatype
from mdvtools.mdvproject import MDVProject


def test_contract_covers_every_legacy_datatype():
    assert set(LEGACY_DATATYPES) == {
        "integer",
        "double",
        "int32",
        "text",
        "text16",
        "unique",
        "multitext",
    }


def test_storage_dtypes_match_what_h5_projects_have_always_used():
    assert NUMPY_STORAGE_DTYPES == {
        "integer": np.float32,
        "double": np.float32,
        "int32": np.int32,
        "text": np.uint8,
        "text16": np.uint16,
        "multitext": np.uint16,
    }


@pytest.mark.parametrize(
    "dtype, expected",
    [
        # numpy
        (np.dtype("int8"), "integer"),
        (np.dtype("int32"), "integer"),
        (np.dtype("int64"), "integer"),
        (np.dtype("uint32"), "integer"),
        (np.dtype("float32"), "double"),
        (np.dtype("float64"), "double"),
        (np.dtype("bool"), "text"),
        (np.dtype("object"), "text"),
        (np.dtype("U5"), "text"),
        (np.dtype("datetime64[ns]"), "double"),
        # pandas extension dtypes, which used to raise KeyError (#480)
        (pd.Int32Dtype(), "integer"),
        (pd.UInt32Dtype(), "integer"),
        (pd.Int64Dtype(), "integer"),
        (pd.Float64Dtype(), "double"),
        (pd.BooleanDtype(), "text"),
        (pd.StringDtype(), "text"),
        (pd.CategoricalDtype(["a"]), "text"),
        (pd.DatetimeTZDtype(tz="UTC"), "double"),
        # polars, classes and instances
        (pl.Int32, "integer"),
        (pl.UInt64(), "integer"),
        (pl.Float32, "double"),
        (pl.Boolean, "text"),
        (pl.String, "text"),
        (pl.Categorical, "text"),
        (pl.Date, "double"),
        (pl.Datetime("us"), "double"),
    ],
)
def test_infer_datatype(dtype, expected):
    assert infer_datatype(dtype) == expected


def test_add_datasource_accepts_pandas_nullable_dtypes(tmp_path):
    """#480: nullable extension dtypes failed with KeyError in get_column_info."""
    df = pd.DataFrame(
        {
            "count": pd.array([1, None, 3], dtype="Int32"),
            "score": pd.array([0.5, None, 1.5], dtype="Float64"),
            "flag": pd.array([True, None, False], dtype="boolean"),
        }
    )
    project = MDVProject(str(tmp_path / "project"), delete_existing=True)

    project.add_datasource("cells", df, add_to_view=None)

    columns = {c["field"]: c for c in project.get_datasource_metadata("cells")["columns"]}
    assert columns["count"]["datatype"] == "integer"
    assert columns["score"]["datatype"] == "double"
    assert columns["flag"]["datatype"] == "text"
    count = project.get_column("cells", "count")
    assert count[0] == 1 and np.isnan(count[1]) and count[2] == 3
