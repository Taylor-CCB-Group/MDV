"""Tests for datetime → days-since-epoch date column ingest."""

from __future__ import annotations

import numpy as np
import pandas as pd
import polars as pl
import pytest

from mdvtools.mdvproject import (
    MDVProject,
    apply_date_column_metadata,
    convert_pandas_datetime_columns,
    map_polars_datetime_to_day_doubles,
)


def test_convert_pandas_datetime_columns_sets_day_doubles():
    df = pd.DataFrame(
        {
            "when": pd.to_datetime(
                ["1970-01-01", "2020-01-01", None],
                utc=True,
            ),
            "label": ["a", "b", "c"],
        }
    )
    out, date_fields = convert_pandas_datetime_columns(df)
    assert date_fields == {"when"}
    assert out["when"].dtype == np.float64
    assert out["when"].iloc[0] == pytest.approx(0.0)
    assert out["when"].iloc[1] == pytest.approx(18262.0)
    assert np.isnan(out["when"].iloc[2])
    assert list(out["label"]) == ["a", "b", "c"]


def test_convert_pandas_object_iso_dates():
    df = pd.DataFrame({"d": ["2024-06-01", "2024-06-15", "2024-07-01"]})
    out, date_fields = convert_pandas_datetime_columns(df)
    assert date_fields == {"d"}
    assert out["d"].iloc[0] < out["d"].iloc[1] < out["d"].iloc[2]
    span = out["d"].iloc[2] - out["d"].iloc[0]
    assert span == pytest.approx(30.0)


def test_apply_date_column_metadata():
    cols = [
        {"field": "when", "datatype": "double"},
        {"field": "x", "datatype": "double"},
    ]
    apply_date_column_metadata(cols, {"when"})
    assert cols[0]["is_date"] is True
    assert cols[0]["date_unit"] == "days"
    assert "is_date" not in cols[1]


def test_add_datasource_datetime_column(tmp_path):
    project = MDVProject(str(tmp_path), delete_existing=True)
    df = pd.DataFrame(
        {
            "when": pd.to_datetime(["2020-01-01", "2020-01-11", "2020-02-01"]),
            "value": [1.0, 2.0, 3.0],
        }
    )
    project.add_datasource("events", df, add_to_view=None)
    meta = project.get_datasource_metadata("events")
    when_col = next(c for c in meta["columns"] if c["field"] == "when")
    assert when_col["datatype"] == "double"
    assert when_col["is_date"] is True
    assert when_col["date_unit"] == "days"
    assert when_col["minMax"][1] - when_col["minMax"][0] == pytest.approx(31.0)


def test_polars_datetime_to_day_doubles():
    s = pl.Series("d", ["1970-01-01", "2020-01-01"]).str.to_date()
    days = map_polars_datetime_to_day_doubles(s)
    assert days[0] == pytest.approx(0.0)
    assert days[1] == pytest.approx(18262.0)


def test_polars_is_date_on_already_numeric_preserves_day_doubles(tmp_path):
    """is_date metadata on Float64 day doubles must not re-run datetime conversion."""
    project = MDVProject(str(tmp_path), delete_existing=True)
    df = pl.DataFrame(
        {
            "when": [0.0, 18262.0, 18300.0],
            "value": [1.0, 2.0, 3.0],
        }
    )
    columns = [
        {
            "field": "when",
            "name": "when",
            "datatype": "double",
            "is_date": True,
            "date_unit": "days",
        },
        {"field": "value", "name": "value", "datatype": "double"},
    ]
    project.add_datasource_polars("events", df, columns=columns, add_to_view=None)
    meta = project.get_datasource_metadata("events")
    when_col = next(c for c in meta["columns"] if c["field"] == "when")
    assert when_col["is_date"] is True
    assert when_col["date_unit"] == "days"
    assert when_col["minMax"][0] == pytest.approx(0.0)
    assert when_col["minMax"][1] == pytest.approx(18300.0)
