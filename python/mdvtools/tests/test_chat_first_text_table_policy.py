from types import SimpleNamespace
from typing import cast

from langchain.prompts import PromptTemplate

from mdvtools.llm import templates
from mdvtools.llm.verification import build_verification_summary
from mdvtools.mdvproject import MDVProject


class _FakeProject:
    def __init__(self):
        self.datasources = [{"name": "cells"}]
        self.dir = "/tmp/fake-project"

    def get_datasource_metadata(self, name):
        if name == "cells":
            return {"name": "cells", "size": 10, "columns": [{"field": "leiden", "name": "Leiden"}]}
        raise KeyError(name)

    def get_view(self, _name):
        return {"initialCharts": {}}


def test_prompt_includes_chat_first_text_table_policy(monkeypatch):
    fake_roles = SimpleNamespace(
        expressions=[],
        obs_datasource="cells",
    )

    monkeypatch.setattr(templates, "infer_datasource_roles", lambda _project: fake_roles)
    monkeypatch.setattr(templates, "create_column_markdown", lambda _cols: "columns-md")

    prompt = templates.get_createproject_prompt_RAG(
        project=cast(MDVProject, _FakeProject()),
        path_to_data="",
        datasource_name="cells",
        final_answer='fields "leiden"\ncharts "Table Plot"',
        question="list predicted cell types by cluster",
    )

    assert "Chat-first textual/table outputs" in prompt
    assert "Precedence vs marker persistence" in prompt
    assert "interpretation" in prompt.lower()
    assert "**narrow** ChatMDV" in prompt or "narrow** ChatMDV" in prompt
    assert "chat_rank_genes_result" in prompt
    assert "decide **whether** a saved-view chart is needed" in prompt
    assert "print/markdown-only" in prompt
    assert "Do NOT create `TextBox` or `TablePlot` by default" in prompt
    assert "Do not add a selection dialog unconditionally." in prompt
    assert "Feature table field compatibility" in prompt
    assert "Field ID" in prompt
    assert "Marker genes and missing columns" in prompt
    assert "Visualization vs analysis consistency" in prompt
    assert "Single source of truth" in prompt
    assert "Marker ranking vs DotPlot" in prompt
    assert "Hybrid Scanpy routing contract" in prompt
    assert "ID-safe merge required for `set_column(...)`" in prompt
    assert "Pre-chart validation gate" in prompt
    assert "those genes" in prompt
    assert "heatmap/dot/bubble/violin" in prompt.lower()
    assert "cluster count/distribution" in prompt.lower()
    assert "cell-type prediction" in prompt.lower()
    assert "build `params` from persisted datasource metadata" in prompt
    assert "skip saving a broken chart" in prompt
    assert "subgroup_key" in prompt
    assert "rna_expr" in prompt.lower()
    assert "Metadata-first chart params" in prompt
    assert "Never** invent" in prompt or "invent a datasource" in prompt


def test_prompt_includes_marker_gene_policy_with_h5ad(monkeypatch):
    fake_roles = SimpleNamespace(
        expressions=[],
        obs_datasource="cells",
    )

    monkeypatch.setattr(templates, "infer_datasource_roles", lambda _project: fake_roles)
    monkeypatch.setattr(templates, "create_column_markdown", lambda _cols: "columns-md")

    prompt = templates.get_createproject_prompt_RAG(
        project=cast(MDVProject, _FakeProject()),
        path_to_data="/project/scanpy_data.h5ad",
        datasource_name="cells",
        final_answer='fields "leiden"\ncharts "Table Plot"',
        question="top 5 marker genes per cluster",
    )

    assert "rank_genes_groups" in prompt
    assert "Marker genes and missing columns" in prompt
    assert "Visualization vs analysis consistency" in prompt
    assert "When user says “use scanpy”" in prompt or 'When user says “use scanpy”' in prompt


def test_createproject_prompt_RAG_formats_as_retrieval_qa_template(monkeypatch):
    """RAG body must not contain single-brace placeholders except {context} (LangChain PromptTemplate)."""
    fake_roles = SimpleNamespace(
        expressions=[],
        obs_datasource="cells",
    )

    monkeypatch.setattr(templates, "infer_datasource_roles", lambda _project: fake_roles)
    monkeypatch.setattr(templates, "create_column_markdown", lambda _cols: "columns-md")

    prompt = templates.get_createproject_prompt_RAG(
        project=cast(MDVProject, _FakeProject()),
        path_to_data="",
        datasource_name="cells",
        final_answer='fields "leiden"\ncharts "Table Plot"',
        question="top markers",
    )

    pt = PromptTemplate(template=prompt, input_variables=["context", "question"])
    out = pt.format(context="CTX", question="Q")
    assert "CTX" in out


def test_verification_uses_response_wording_not_view_wording():
    summary = build_verification_summary(
        project=_FakeProject(),
        final_code='datasource_name = "cells"\nproject.set_view("v", {"initialCharts": {}})\n',
        view_name="v",
    )

    assert "not used in this response" in summary
    assert "### Columns / outputs (from code)" in summary
    assert "### Chart types (if any, from code)" in summary
