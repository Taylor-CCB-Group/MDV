from mdvtools.llm.code_preflight import validate_generated_code_preflight
from mdvtools.llm.preflight_flow import preflight_with_single_retry


def test_preflight_detects_missing_import_for_chart_class():
    code = """
from mdvtools.charts.box_plot import BoxPlot

def main():
    plot = MultilineChart(title="x", params=["a", "b"])
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is False
    msgs = " | ".join(i.message for i in res.issues)
    assert "Unknown chart class `MultilineChart`" in msgs


def test_preflight_rejects_unknown_chart_class():
    code = """
def main():
    c = RadarChart(title="x", params=["a"])
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is False
    assert any(i.code == "unknown_chart_class" for i in res.issues)


def test_preflight_allows_valid_multi_chart_script():
    code = """
from mdvtools.charts.box_plot import BoxPlot
from mdvtools.charts.violin_plot import ViolinPlot
from mdvtools.charts.multi_line_plot import MultiLinePlot

def main():
    b = BoxPlot(title="b", params=["leiden", "n_genes"], size=[1, 1], position=[0, 0])
    v = ViolinPlot(title="v", params=["leiden", "n_genes"], size=[1, 1], position=[0, 0])
    m = MultiLinePlot(title="m", params=["leiden", "n_genes"], size=[1, 1], position=[0, 0])
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is True
    assert res.issues == []
    assert "BoxPlot" in res.chart_classes
    assert "ViolinPlot" in res.chart_classes
    assert "MultiLinePlot" in res.chart_classes


def test_preflight_flow_retries_once_and_returns_metadata():
    bad = """
def main():
    c = MultilineChart(title="x", params=["a", "b"])
"""
    good = """
from mdvtools.charts.multi_line_plot import MultiLinePlot
def main():
    c = MultiLinePlot(title="x", params=["a", "b"])
"""
    calls: list[str] = []

    def regenerate(issue_text: str) -> str:
        calls.append(issue_text)
        return good

    final_code, meta = preflight_with_single_retry(
        initial_code=bad,
        regenerate_once=regenerate,
        log=lambda *_: None,
    )
    assert final_code == good
    assert meta["preflight_retried"] is True
    assert meta["preflight_ok"] is True
    assert len(calls) == 1


def test_preflight_detects_invalid_column_for_datasource():
    code = """
from mdvtools.mdvproject import MDVProject

def main():
    datasource_name = "cells"
    columns = ["bucket", "name"]
    project = MDVProject("/tmp/project", delete_existing=False)
    df = project.get_datasource_as_dataframe(datasource_name, columns=columns)
"""
    res = validate_generated_code_preflight(
        code,
        datasource_fields={
            "cells": {"bucket", "cell_id", "major"},
            "genes": {"name", "gene_id"},
        },
    )
    assert res.ok is False
    assert any(i.code == "unknown_datasource_column" for i in res.issues)
    assert any("Datasource `cells` does not contain column(s) ['name']" in i.message for i in res.issues)


def test_preflight_flow_retries_on_invalid_datasource_columns():
    bad = """
from mdvtools.mdvproject import MDVProject

def main():
    datasource_name = "cells"
    columns = ["bucket", "name"]
    project = MDVProject("/tmp/project", delete_existing=False)
    df = project.get_datasource_as_dataframe(datasource_name, columns=columns)
"""
    good = """
from mdvtools.mdvproject import MDVProject

def main():
    datasource_name = "cells"
    columns = ["bucket", "cell_id"]
    project = MDVProject("/tmp/project", delete_existing=False)
    df = project.get_datasource_as_dataframe(datasource_name, columns=columns)
"""
    calls: list[str] = []

    def regenerate(issue_text: str) -> str:
        calls.append(issue_text)
        return good

    final_code, meta = preflight_with_single_retry(
        initial_code=bad,
        regenerate_once=regenerate,
        log=lambda *_: None,
        datasource_fields={"cells": {"bucket", "cell_id"}},
    )
    assert final_code == good
    assert meta["preflight_retried"] is True
    assert meta["preflight_ok"] is True
    assert len(calls) == 1

