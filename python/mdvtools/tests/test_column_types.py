"""The column-type contract (ADR-0024) as the Python package reads it."""

import numpy as np
import pandas as pd
import polars as pl
import pytest

from mdvtools.column_types import (
    LEGACY_DATATYPES,
    NUMPY_STORAGE_DTYPES,
    infer_datatype,
    integer_precision_warning,
)
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
        # numpy integers that fit int32 are stored exactly, as int32 (#537)
        (np.dtype("int8"), "int32"),
        (np.dtype("int16"), "int32"),
        (np.dtype("int32"), "int32"),
        (np.dtype("uint8"), "int32"),
        (np.dtype("uint16"), "int32"),
        # wider ones go to float32-backed integer
        (np.dtype("int64"), "integer"),
        (np.dtype("uint32"), "integer"),
        (np.dtype("uint64"), "integer"),
        (np.dtype("float32"), "double"),
        (np.dtype("float64"), "double"),
        (np.dtype("bool"), "text"),
        (np.dtype("object"), "text"),
        (np.dtype("U5"), "text"),
        (np.dtype("datetime64[ns]"), "double"),
        # pandas nullable dtypes, which used to raise KeyError (#480); int32 storage has no
        # missing value, so nullable integers go to float32-backed integer
        (pd.Int8Dtype(), "integer"),
        (pd.Int32Dtype(), "integer"),
        (pd.UInt32Dtype(), "integer"),
        (pd.Int64Dtype(), "integer"),
        (pd.Float64Dtype(), "double"),
        (pd.BooleanDtype(), "text"),
        (pd.StringDtype(), "text"),
        (pd.CategoricalDtype(["a"]), "text"),
        (pd.DatetimeTZDtype(tz="UTC"), "double"),
        # polars, classes and instances; every polars column can hold nulls
        (pl.Int8, "integer"),
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


@pytest.mark.parametrize("original_dtype", ["int32", "Int32", "int64", "Int64", "uint32", "UInt32", "UInt64"])
def test_integer_precision_warning_for_wide_integers(original_dtype):
    warning = integer_precision_warning(original_dtype, np.array([1, -16_777_217, np.nan]))
    assert warning is not None
    assert warning.startswith(f"{original_dtype} -> float32 (integer): 1 value(s)")


def test_no_integer_precision_warning_when_values_fit_or_type_is_narrow():
    assert integer_precision_warning("int64", np.array([1, 16_777_216])) is None
    assert integer_precision_warning("int16", np.array([1, 2])) is None
    assert integer_precision_warning("float64", np.array([1e9])) is None


def test_add_datasource_stores_int32_columns_exactly(tmp_path):
    """A float32 cannot hold 16_777_217; an int32 column can."""
    df = pd.DataFrame({"id": np.array([1, 16_777_217], dtype="int32")})
    project = MDVProject(str(tmp_path / "project"), delete_existing=True)

    project.add_datasource("cells", df, add_to_view=None)

    column = project.get_datasource_metadata("cells")["columns"][0]
    assert column["datatype"] == "int32"
    assert list(project.get_column("cells", "id")) == [1, 16_777_217]


@pytest.mark.parametrize("dtype", ["UInt32", "int64"])
def test_add_datasource_warns_when_wide_integers_lose_precision(tmp_path, dtype):
    df = pd.DataFrame({"id": pd.array([1, 16_777_217], dtype=dtype)})
    project = MDVProject(str(tmp_path / "project"), delete_existing=True)

    project.add_datasource("cells", df, add_to_view=None)

    column = project.get_datasource_metadata("cells")["columns"][0]
    assert column["datatype"] == "integer"
    assert any("precision loss possible" in w for w in column.get("storage_warnings", []))


def test_add_datasource_polars_warns_when_wide_integers_lose_precision(tmp_path):
    df = pl.DataFrame({"id": pl.Series([1, 16_777_217], dtype=pl.Int64)})
    project = MDVProject(str(tmp_path / "project"), delete_existing=True)

    project.add_datasource_polars("cells", df, add_to_view=None)

    column = project.get_datasource_metadata("cells")["columns"][0]
    assert column["datatype"] == "integer"
    assert any("precision loss possible" in w for w in column.get("storage_warnings", []))
