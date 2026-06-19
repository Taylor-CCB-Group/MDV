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


def test_preflight_rejects_unknown_histogram_alias():
    code = """
def main():
    h = Histogram(title="h", params=["Age"])
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is False
    msgs = " | ".join(i.message for i in res.issues)
    assert "Unknown chart class `Histogram`" in msgs
    assert "did you mean `HistogramPlot`?" in msgs


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


def test_preflight_allows_histogram_plot_class():
    code = """
from mdvtools.charts.histogram_plot import HistogramPlot

def main():
    h = HistogramPlot(
        title="h",
        param="Age",
        bin_number=10,
        display_min=0,
        display_max=100,
        size=[1, 1],
        position=[0, 0],
    )
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is True
    assert res.issues == []
    assert "HistogramPlot" in res.chart_classes


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


def test_preflight_hints_owner_datasource_for_missing_column():
    code = """
from mdvtools.mdvproject import MDVProject

def main():
    project = MDVProject("/tmp/project", delete_existing=False)
    df = project.get_datasource_as_dataframe("qc_channels", columns=["channel_name", "cv_pct"])
"""
    res = validate_generated_code_preflight(
        code,
        datasource_fields={
            "qc_channels": {"channel_name", "run_id"},
            "qc_field_uniformity": {"channel_name", "cv_pct", "uniformity_score"},
        },
    )
    assert res.ok is False
    msgs = " | ".join(i.message for i in res.issues)
    assert "cv_pct" in msgs
    assert "qc_field_uniformity" in msgs


def test_preflight_rejects_set_row_indices_on_scatter_plot():
    code = """
from mdvtools.charts.scatter_plot import ScatterPlot

def main():
    plot = ScatterPlot(title="t", params=["x", "y"], size=[1, 1], position=[0, 0])
    plot.set_row_indices([0, 1, 2])
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is False
    assert any(i.code == "forbidden_chart_method" for i in res.issues)
    assert any("set_row_indices" in i.message for i in res.issues)


def test_preflight_rejects_set_background_filter_even_without_chart_binding():
    code = """
from mdvtools.charts.scatter_plot import ScatterPlot

def main():
    plot = ScatterPlot(title="t", params=["x", "y"], size=[1, 1], position=[0, 0])
    plot.set_background_filter({"column": "major", "values": ["B"]})
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is False
    assert any(i.code == "forbidden_chart_method" for i in res.issues)
    assert any("set_background_filter" in i.message for i in res.issues)


def test_preflight_allows_scatter_plot_setters():
    code = """
from mdvtools.charts.scatter_plot import ScatterPlot

def main():
    plot = ScatterPlot(title="t", params=["x", "y"], size=[1, 1], position=[0, 0])
    plot.set_color_by("expr|GENE(expr)|0")
    plot.set_axis_properties("x", {"label": "UMAP1"})
    plot.set_brush("poly")
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is True
    assert res.issues == []


def test_preflight_rejects_unknown_method_on_tracked_chart():
    code = """
from mdvtools.charts.scatter_plot import ScatterPlot

def main():
    plot = ScatterPlot(title="t", params=["x", "y"], size=[1, 1], position=[0, 0])
    plot.subset_by_category("major", ["B"])
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is False
    assert any(i.code == "unknown_chart_method" for i in res.issues)


def test_preflight_flow_retries_on_forbidden_chart_method():
    bad = """
from mdvtools.charts.scatter_plot import ScatterPlot

def main():
    plot = ScatterPlot(title="t", params=["x", "y"], size=[1, 1], position=[0, 0])
    plot.set_row_indices([0, 1])
"""
    good = """
from mdvtools.charts.scatter_plot import ScatterPlot

def main():
    plot = ScatterPlot(title="t", params=["x", "y"], size=[1, 1], position=[0, 0])
    plot.set_color_by("expr|GENE(expr)|0")
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
    assert "set_row_indices" in calls[0]


def test_preflight_rejects_col_dtype_on_metadata_columns():
    code = """
from mdvtools.mdvproject import MDVProject

def main():
    project = MDVProject("/tmp/project", delete_existing=False)
    obs_metadata = project.get_datasource_metadata("cells")
    categorical_fields = [col['name'] for col in obs_metadata['columns'] if col['dtype'] == 'text']
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is False
    assert any(i.code == "invalid_metadata_column_key" for i in res.issues)
    assert any("datatype" in i.message for i in res.issues)


def test_preflight_rejects_col_get_dtype():
    code = """
def main():
    col = {"field": "x", "datatype": "text"}
    if col.get('dtype') == 'text':
        pass
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is False
    assert any(i.code == "invalid_metadata_column_key" for i in res.issues)


def test_preflight_allows_col_get_datatype():
    code = """
from mdvtools.llm.datasource_roles import categorical_field_ids_from_metadata
from mdvtools.mdvproject import MDVProject

def main():
    project = MDVProject("/tmp/project", delete_existing=False)
    obs_metadata = project.get_datasource_metadata("cells")
    categorical_fields = categorical_field_ids_from_metadata(obs_metadata)
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is True


def test_preflight_flow_retries_on_col_dtype():
    bad = """
def main():
    obs_metadata = {"columns": [{"field": "major", "datatype": "text", "name": "Major"}]}
    cats = [c["field"] for c in obs_metadata["columns"] if c["dtype"] == "text"]
"""
    good = """
from mdvtools.llm.datasource_roles import categorical_field_ids_from_metadata

def main():
    obs_metadata = {"columns": [{"field": "major", "datatype": "text", "name": "Major"}]}
    cats = categorical_field_ids_from_metadata(obs_metadata)
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
    assert meta["preflight_ok"] is True
    assert "dtype" in calls[0]


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


def test_preflight_ignores_wrapper_literal_with_trailing_garbage():
    code = '''
def main():
    _ = "bad|MS4A1(bad)|3__junk"
'''
    res = validate_generated_code_preflight(
        code, allowed_wrapper_subgroup_keys={"gs"}
    )
    assert not any(i.code == "invalid_wrapper_subgroup" for i in res.issues)

