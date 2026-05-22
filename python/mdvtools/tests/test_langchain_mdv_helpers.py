from mdvtools.llm.langchain_mdv import (
    _is_text_table_intent,
    _is_text_table_only_code,
)


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
