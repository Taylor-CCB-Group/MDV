from mdvtools.llm.column_field_resolve import (
    normalize_view_chart_params,
    prune_view_charts_with_invalid_params,
)


class FakeProjectCells:
    """cells has leiden only; no gene_ids."""

    def get_datasource_metadata(self, name: str):
        assert name == "cells"
        return {
            "name": "cells",
            "columns": [
                {"name": "leiden", "field": "leiden", "datatype": "text"},
            ],
        }


def test_prune_drops_table_chart_with_unknown_param_tokens():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "markers",
                    "type": "table_chart",
                    "param": ["leiden", "gene_ids", "score"],
                }
            ]
        }
    }
    proj = FakeProjectCells()
    pruned, n = prune_view_charts_with_invalid_params(view, proj)
    assert n == 1
    assert pruned["initialCharts"]["cells"] == []


def test_prune_keeps_wrapper_params():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "x",
                    "type": "dot_plot",
                    "param": ["rna_expr|CD3D(rna_expr)|12", "leiden"],
                }
            ]
        }
    }
    proj = FakeProjectCells()
    pruned, n = prune_view_charts_with_invalid_params(view, proj)
    assert n == 0
    assert len(pruned["initialCharts"]["cells"]) == 1


def test_normalize_then_prune_removes_bad_marker_chart():
    view = {
        "initialCharts": {
            "cells": [
                {
                    "title": "markers",
                    "type": "table_chart",
                    "param": ["leiden", "gene_ids"],
                }
            ]
        }
    }
    proj = FakeProjectCells()
    normalized = normalize_view_chart_params(view, proj)
    pruned, n = prune_view_charts_with_invalid_params(normalized, proj)
    assert n == 1
    assert pruned["initialCharts"]["cells"] == []
