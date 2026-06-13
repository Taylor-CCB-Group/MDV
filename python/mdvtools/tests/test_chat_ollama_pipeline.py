"""Tests for Ollama-specific ChatMDV pipeline guards and compact RAG prompts."""

from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from mdvtools.llm import templates
from mdvtools.llm.langchain_mdv import (
    ProjectChat,
    _agent_fields_charts_output_valid,
    _fallback_agent_plan,
    _infer_chart_types_from_question,
    _invoke_rag_with_empty_retry,
    _raise_if_no_extracted_code,
    _resolve_agent_plan,
)
from mdvtools.llm.llm_providers import ModelSpec


class _FakeProject:
    def __init__(self):
        self.datasources = [{"name": "cells"}]
        self.dir = "/tmp/fake-project"

    def get_datasource_metadata(self, name: str):
        if name == "cells":
            return {"name": "cells", "size": 10, "columns": [{"field": "mean_counts", "name": "mean_counts"}]}
        raise KeyError(name)

    def get_datasource_names(self) -> list[str]:
        return [str(ds["name"]) for ds in self.datasources]


def test_agent_fields_charts_output_valid_accepts_plan():
    assert _agent_fields_charts_output_valid('fields "mean_counts"\ncharts "Histogram"') is True


def test_agent_fields_charts_output_valid_rejects_empty():
    assert _agent_fields_charts_output_valid("") is False
    assert _agent_fields_charts_output_valid("   ") is False
    assert _agent_fields_charts_output_valid("no plan here") is False
    assert _agent_fields_charts_output_valid("Agent stopped due to iteration limit or time limit.") is False


def test_resolve_agent_plan_from_react_intermediate_steps():
    action = SimpleNamespace(
        log='fields "X_umap_1", "X_umap_2", "leiden"\ncharts "Scatter plot (2D)"',
        tool_input="Invalid Format: Missing 'Action:' after 'Thought:'",
    )
    response = {
        "output": "Agent stopped due to iteration limit or time limit.",
        "intermediate_steps": [(action, "obs")],
    }
    plan = _resolve_agent_plan(response)
    assert 'fields "X_umap_1"' in plan
    assert 'charts "Scatter plot (2D)"' in plan
    assert _agent_fields_charts_output_valid(plan)


def test_fallback_agent_plan_from_question():
    plan = _fallback_agent_plan("Show a UMAP plot colored by Leiden clusters.")
    assert _agent_fields_charts_output_valid(plan)
    assert "Scatter plot (2D)" in plan
    assert "infer from Project Data Context" in plan


def test_infer_chart_types_from_question_heatmap():
    assert "Heatmap" in _infer_chart_types_from_question(
        "Compare NKG7 and GNLY across clusters using a heatmap."
    )


def test_create_custom_pandas_agent_uses_functions_first_for_ollama():
    chat = ProjectChat.__new__(ProjectChat)
    captured: dict = {}

    def fake_build(*args, **kwargs):
        captured["use_react"] = kwargs.get("use_react")
        mock_executor = MagicMock()
        mock_executor.invoke.return_value = {"output": 'fields "x"\ncharts "Histogram"'}
        return mock_executor

    with patch.object(ProjectChat, "_build_agent_executor", side_effect=fake_build):
        with patch("mdvtools.llm.langchain_mdv.LLMChain") as mock_chain_cls:
            mock_chain_cls.return_value.run.return_value = "question"
            agent = chat.create_custom_pandas_agent(
                llm=MagicMock(),
                dfs={"df1": MagicMock(), "df2": MagicMock()},
                prompt_data="",
                memory=MagicMock(),
                provider="ollama",
            )
            agent("histogram of mean_counts")

    assert captured.get("use_react") is False


def test_create_custom_pandas_agent_falls_back_to_react_when_ollama_plan_empty():
    chat = ProjectChat.__new__(ProjectChat)
    use_react_calls: list[bool] = []

    def fake_build(*args, **kwargs):
        use_react = kwargs.get("use_react")
        use_react_calls.append(use_react is True)
        mock_executor = MagicMock()
        if use_react:
            action = SimpleNamespace(
                log='fields "leiden"\ncharts "Histogram"',
                tool_input="",
            )
            mock_executor.invoke.return_value = {
                "output": "Agent stopped due to iteration limit or time limit.",
                "intermediate_steps": [(action, "obs")],
            }
        else:
            mock_executor.invoke.return_value = {"output": ""}
        return mock_executor

    with patch.object(ProjectChat, "_build_agent_executor", side_effect=fake_build):
        with patch("mdvtools.llm.langchain_mdv.LLMChain") as mock_chain_cls:
            mock_chain_cls.return_value.run.return_value = "question"
            agent = chat.create_custom_pandas_agent(
                llm=MagicMock(),
                dfs={"df1": MagicMock(), "df2": MagicMock()},
                prompt_data="",
                memory=MagicMock(),
                provider="ollama",
            )
            response = agent("histogram")

    assert use_react_calls == [False, True]
    assert _agent_fields_charts_output_valid(_resolve_agent_plan(response))


def test_create_custom_pandas_agent_uses_functions_for_openai():
    chat = ProjectChat.__new__(ProjectChat)
    captured: dict = {}

    def fake_build(*args, **kwargs):
        captured["use_react"] = kwargs.get("use_react")
        mock_executor = MagicMock()
        mock_executor.invoke.return_value = {"output": 'fields "x"\ncharts "Histogram"'}
        return mock_executor

    with patch.object(ProjectChat, "_build_agent_executor", side_effect=fake_build):
        with patch("mdvtools.llm.langchain_mdv.LLMChain") as mock_chain_cls:
            mock_chain_cls.return_value.run.return_value = "question"
            agent = chat.create_custom_pandas_agent(
                llm=MagicMock(),
                dfs={"df1": MagicMock(), "df2": MagicMock()},
                prompt_data="",
                memory=MagicMock(),
                provider="openai",
            )
            agent("histogram of mean_counts")

    assert captured.get("use_react") is False


def test_invoke_rag_with_empty_retry_raises_after_second_empty():
    chat_model = ModelSpec(
        id="ollama:chat:qwen3-coder:latest",
        label="Qwen3 Coder",
        provider="ollama",
        model="qwen3-coder:latest",
        kind="chat",
        available=True,
    )
    qa_chain = MagicMock()
    qa_chain.invoke.side_effect = [
        {"result": "", "source_documents": []},
        {"result": "   ", "source_documents": []},
    ]

    with pytest.raises(ValueError, match="Code generation returned empty output"):
        _invoke_rag_with_empty_retry(qa_chain, "User question: hist", chat_model, None)

    assert qa_chain.invoke.call_count == 2


def test_invoke_rag_with_empty_retry_succeeds_on_retry():
    chat_model = ModelSpec(
        id="ollama:chat:qwen3-coder:latest",
        label="Qwen3 Coder",
        provider="ollama",
        model="qwen3-coder:latest",
        kind="chat",
        available=True,
    )
    qa_chain = MagicMock()
    qa_chain.invoke.side_effect = [
        {"result": "", "source_documents": []},
        {"result": "```python\nprint('ok')\n```", "source_documents": []},
    ]

    output_qa, result = _invoke_rag_with_empty_retry(
        qa_chain, "User question: hist", chat_model, None
    )
    assert "print('ok')" in result
    assert output_qa["result"].startswith("```python")


def test_raise_if_no_extracted_code_rejects_warning_stub():
    stub = "# WARNING:::: No code captured from the response when calling prepare_code().\n"
    with pytest.raises(ValueError, match="did not produce runnable Python"):
        _raise_if_no_extracted_code("", stub)


def test_raise_if_no_extracted_code_rejects_empty_extraction():
    with pytest.raises(ValueError, match="did not produce runnable Python"):
        _raise_if_no_extracted_code("no fenced code here", "def main(): pass\n")


def test_compact_rag_prompt_is_shorter_than_full(monkeypatch):
    fake_roles = SimpleNamespace(
        expressions=[],
        obs_datasource="cells",
    )
    monkeypatch.setattr(templates, "infer_datasource_roles", lambda _project: fake_roles)
    monkeypatch.setattr(templates, "create_column_markdown", lambda _cols: "columns-md")

    project = _FakeProject()
    kwargs = dict(
        project=project,
        path_to_data="",
        datasource_name="cells",
        final_answer='fields "mean_counts"\ncharts "Histogram"',
        question="Which genes have highest mean_counts?",
    )
    full = templates.get_createproject_prompt_RAG(**kwargs, compact=False)
    compact = templates.get_createproject_prompt_RAG(**kwargs, compact=True)

    assert len(compact) < len(full)
    assert "columns-md" in compact
    assert "Return only one fenced ```python code block" in compact
    assert "Do NOT call `project.add_datasource(...)`" in compact
    assert "Marker ranking vs DotPlot" not in compact
    assert "Visualization vs analysis consistency" not in compact
