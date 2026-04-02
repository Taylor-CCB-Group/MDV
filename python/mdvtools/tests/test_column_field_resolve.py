"""Tests for column_field_resolve (Chat MDV name -> field id)."""

from mdvtools.llm.column_field_resolve import (
    build_name_to_field_map,
    field_set_from_columns,
    normalize_view_chart_params,
    resolve_param_string,
)


def test_build_name_to_field_map_unique():
    cols = [
        {"name": "Foo Bar", "field": "Foo_Bar", "datatype": "double"},
        {"name": "Other", "field": "other_id", "datatype": "text"},
    ]
    m = build_name_to_field_map(cols)
    assert m["Foo Bar"] == "Foo_Bar"
    assert m["Other"] == "other_id"


def test_build_name_to_field_map_ambiguous_omits():
    cols = [
        {"name": "Dup", "field": "a", "datatype": "text"},
        {"name": "Dup", "field": "b", "datatype": "text"},
    ]
    m = build_name_to_field_map(cols)
    assert "Dup" not in m


def test_field_set_from_columns():
    cols = [
        {"name": "X", "field": "fx", "datatype": "double"},
        {"name": "Y", "datatype": "text"},
    ]
    s = field_set_from_columns(cols)
    assert "fx" in s
    assert "Y" in s


def test_resolve_param_string_field_exact():
    fs = {"my_field", "other"}
    nm = {"Display": "my_field"}
    assert resolve_param_string("my_field", fs, nm) == "my_field"


def test_resolve_param_string_name_to_field():
    fs = {"my_field"}
    nm = {"Display Name": "my_field"}
    assert resolve_param_string("Display Name", fs, nm) == "my_field"


def test_resolve_param_string_gs_passthrough():
    fs = set()
    nm = {}
    g = "rna_expr|TP53(rna_expr)|12"
    assert resolve_param_string(g, fs, nm) == g


def test_resolve_param_string_unknown_unchanged():
    fs = {"a"}
    nm = {}
    assert resolve_param_string("not_in_map", fs, nm) == "not_in_map"


def test_normalize_view_chart_params():
    class FakeProject:
        def get_datasource_metadata(self, name):
            assert name == "cells"
            return {
                "name": "cells",
                "columns": [
                    {
                        "name": "Presence of Inflammation",
                        "field": "Presence_of_Inflammation",
                        "datatype": "text",
                    },
                    {
                        "name": "Inflammation Score",
                        "field": "Inflammation_Score",
                        "datatype": "double",
                    },
                ],
            }

    view = {
        "initialCharts": {
            "cells": [
                {
                    "type": "box_plot",
                    "title": "t",
                    "param": ["Presence of Inflammation", "Inflammation Score"],
                    "size": [100, 100],
                    "position": [0, 0],
                    "id": "abc",
                }
            ]
        }
    }
    out = normalize_view_chart_params(view, FakeProject())
    params = out["initialCharts"]["cells"][0]["param"]
    assert params == ["Presence_of_Inflammation", "Inflammation_Score"]
    # original unchanged (deep copy)
    assert view["initialCharts"]["cells"][0]["param"] == [
        "Presence of Inflammation",
        "Inflammation Score",
    ]
