from typing import cast

from mdvtools.llm.langchain_mdv import (
    ProjectChat,
    _is_text_table_intent,
    _is_text_table_only_code,
)
from mdvtools.mdvproject import MDVProject


def test_is_text_table_only_code_print_only():
    assert _is_text_table_only_code("def main():\n    print('ok')\n") is True


def test_is_text_table_only_code_allows_textbox_tableplot():
    code = "from mdvtools.charts.text_box_plot import TextBox\nTextBox(title='t', params=[])"
    assert _is_text_table_only_code(code) is True


def test_is_text_table_only_code_rejects_scatter():
    code = "from mdvtools.charts.scatter_plot import ScatterPlot\nScatterPlot(title='s', params=['x'])"
    assert _is_text_table_only_code(code) is False


def test_text_table_intent_skips_retry_when_output_already_text_table_only():
    question = "list predicted cell types by cluster"
    chart_code = "from mdvtools.charts.scatter_plot import ScatterPlot\nScatterPlot(title='s', params=['x'])"
    text_code = "def main():\n    print(df.head())\n"
    assert _is_text_table_intent(question) is True
    assert _is_text_table_intent(question) and not _is_text_table_only_code(chart_code)
    assert not (_is_text_table_intent(question) and not _is_text_table_only_code(text_code))


def test_project_chat_initialization_fails_without_openai_key(monkeypatch):
    class FakeProject:
        datasources = [{"name": "cells"}]

        def get_datasource_as_dataframe(self, _name):
            raise AssertionError("datasource loading should not run without an API key")

    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.setattr("mdvtools.llm.langchain_mdv.load_dotenv", lambda: False)

    def _no_providers() -> None:
        raise ValueError(
            "No LLM providers are available. Set OPENAI_API_KEY or start Ollama "
            "at http://localhost:11434."
        )

    monkeypatch.setattr(
        "mdvtools.llm.langchain_mdv.ensure_any_provider_available",
        _no_providers,
    )

    chat = ProjectChat(cast(MDVProject, FakeProject()))

    assert chat.init_error is True
    assert chat.error_message is not None
    assert "No LLM providers are available" in chat.error_message
