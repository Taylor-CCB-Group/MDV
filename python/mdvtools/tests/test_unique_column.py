"""Tests for the 'unique' column datatype in MDVProject.

Covers round-trip storage, stringLength validation/auto-expansion,
numeric rejection, numpy array inputs, and edge cases.
"""

import os
import tempfile
import shutil
import pytest
import numpy
import pandas

from mdvtools.mdvproject import MDVProject


@pytest.fixture()
def project(tmp_path):
    """Create a throwaway MDVProject with a small datasource."""
    p = MDVProject(str(tmp_path / "proj"), delete_existing=True)
    df = pandas.DataFrame({"id": range(5)})
    p.add_datasource("ds", df)
    return p


# ── round-trip via add_datasource (auto-detection path) ──────────────

def test_unique_column_round_trip_via_set_column(project):
    """set_column with many unique strings should produce a readable unique column."""
    values = [f"item_{i}" for i in range(5)]
    project.set_column("ds", {"name": "labels", "datatype": "unique"}, values)

    meta = project.get_column_metadata("ds", "labels")
    assert meta["datatype"] == "unique"
    assert "stringLength" in meta and meta["stringLength"] > 0

    readback = project.get_column("ds", "labels")
    assert readback == values


# ── set_column_with_raw_data: valid strings ──────────────────────────

def test_set_column_with_raw_data_valid_strings(project):
    """Basic happy path: list of strings with correct stringLength."""
    values = ["alpha", "beta", "gamma", "delta", "epsilon"]
    max_len = max(len(s.encode("utf-8")) for s in values)
    col = {
        "field": "words",
        "name": "words",
        "datatype": "unique",
        "stringLength": max_len,
    }
    project.set_column_with_raw_data("ds", col, values)

    readback = project.get_column("ds", "words")
    assert readback == values


# ── stringLength auto-expansion ──────────────────────────────────────

def test_string_length_auto_expansion(project):
    """If a value exceeds the declared stringLength the metadata must grow."""
    values = ["hi", "hello_world_this_is_long"]
    col = {
        "field": "stretch",
        "name": "stretch",
        "datatype": "unique",
        "stringLength": 2,
    }
    project.set_column_with_raw_data("ds", col, values)

    assert col["stringLength"] >= len("hello_world_this_is_long".encode("utf-8"))
    readback = project.get_column("ds", "stretch")
    assert readback == values


# ── longer ASCII strings ─────────────────────────────────────────────

def test_varying_length_strings_round_trip(project):
    """Strings of different lengths must all survive storage and retrieval."""
    values = ["a", "ab", "abc", "abcdefghij", "z"]
    max_byte_len = max(len(s.encode("utf-8")) for s in values)
    col = {
        "field": "varied",
        "name": "varied",
        "datatype": "unique",
        "stringLength": max_byte_len,
    }
    project.set_column_with_raw_data("ds", col, values)
    readback = project.get_column("ds", "varied")
    assert readback == values


# ── set_column_with_raw_data: no conversion (caller sends strings) ───

def test_set_column_with_raw_data_rejects_int(project):
    """Raw API expects strings; use set_column() for int-to-string (e.g. CSV cell_id)."""
    col = {"field": "bad", "name": "bad", "datatype": "unique", "stringLength": 10}
    with pytest.raises(ValueError, match="expects array of strings|Use set_column"):
        project.set_column_with_raw_data("ds", col, [1, 2, 3])


def test_set_column_with_raw_data_rejects_float(project):
    col = {"field": "bad", "name": "bad", "datatype": "unique", "stringLength": 10}
    with pytest.raises(ValueError, match="expects array of strings|Use set_column"):
        project.set_column_with_raw_data("ds", col, [1.0, 2.5])


def test_set_column_with_raw_data_rejects_numpy_integer_array(project):
    col = {"field": "bad", "name": "bad", "datatype": "unique", "stringLength": 10}
    with pytest.raises(ValueError, match="expects array of strings|Use set_column"):
        project.set_column_with_raw_data("ds", col, numpy.array([1, 2, 3]))


def test_set_column_with_raw_data_rejects_mixed_string_then_numeric(project):
    col = {"field": "bad", "name": "bad", "datatype": "unique", "stringLength": 10}
    with pytest.raises(ValueError, match="expects array of strings|Use set_column"):
        project.set_column_with_raw_data("ds", col, ["ok", 42, "also_ok"])


# ── set_column: int/numeric converted to string for unique ─────────────

def test_set_column_unique_accepts_int_converts_to_string(project):
    """User adding column with datatype 'unique' via set_column (e.g. CSV cell_id) gets string storage."""
    project.set_column(
        "ds",
        {"name": "cell_id", "field": "cell_id", "datatype": "unique"},
        [1, 2, 3, 4, 5],
    )
    readback = project.get_column("ds", "cell_id")
    assert readback == ["1", "2", "3", "4", "5"]


def test_set_column_unique_accepts_mixed_str_and_int(project):
    """set_column with unique and mixed str/int converts ints to strings."""
    project.set_column(
        "ds",
        {"name": "id", "field": "id", "datatype": "unique"},
        ["a", 2, "c", 0],
    )
    readback = project.get_column("ds", "id")
    assert readback == ["a", "2", "c", "0"]


# ── numpy string array input (exercises the ndarray truthiness fix) ──

def test_numpy_string_array(project):
    """numpy.ndarray of strings previously blew up on `if raw_data`.

    h5py cannot write numpy Unicode arrays to fixed-length UTF-8 datasets,
    so we convert to list here.  The important thing is that the validation
    code (len() instead of bool()) does not raise ValueError on ndarray.
    """
    np_values = numpy.array(["foo", "bar", "baz", "qux", "quux"])
    values = list(np_values)
    max_len = max(len(s.encode("utf-8")) for s in values)
    col = {
        "field": "np_str",
        "name": "np_str",
        "datatype": "unique",
        "stringLength": max_len,
    }
    project.set_column_with_raw_data("ds", col, values)
    readback = project.get_column("ds", "np_str")
    assert readback == values


# ── None handling ────────────────────────────────────────────────────

def test_none_values_in_raw_data(project):
    """None should be treated as zero-length, not crash."""
    values = ["a", None, "b", None, "c"]
    col = {
        "field": "nones",
        "name": "nones",
        "datatype": "unique",
        "stringLength": 5,
    }
    project.set_column_with_raw_data("ds", col, values)
    readback = project.get_column("ds", "nones")
    assert readback == ["a", "", "b", "", "c"]


# ── missing / invalid stringLength ───────────────────────────────────

def test_missing_string_length_raises(project):
    col = {"field": "bad", "name": "bad", "datatype": "unique"}
    with pytest.raises(ValueError, match="requires 'stringLength'"):
        project.set_column_with_raw_data("ds", col, ["x"])


def test_zero_string_length_raises(project):
    col = {"field": "bad", "name": "bad", "datatype": "unique", "stringLength": 0}
    with pytest.raises(ValueError, match="requires 'stringLength'"):
        project.set_column_with_raw_data("ds", col, ["x"])


def test_negative_string_length_raises(project):
    col = {"field": "bad", "name": "bad", "datatype": "unique", "stringLength": -5}
    with pytest.raises(ValueError, match="requires 'stringLength'"):
        project.set_column_with_raw_data("ds", col, ["x"])


# ── unknown datatype ─────────────────────────────────────────────────

def test_unknown_datatype_raises(project):
    col = {"field": "bad", "name": "bad", "datatype": "invented_type"}
    with pytest.raises(ValueError, match="Unknown datatype"):
        project.set_column_with_raw_data("ds", col, ["x"])
