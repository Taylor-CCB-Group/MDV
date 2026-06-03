from mdvtools.llm.langchain_mdv import _build_rag_retrieval_query


def test_build_rag_retrieval_query_combines_question_and_charts():
    q = "What is the relationship between total_counts and pct_counts_mt?"
    charts = '"scatter plot (2D)", "table plot"'
    out = _build_rag_retrieval_query(q, charts)
    assert "User question:" in out
    assert q in out
    assert "Suggested chart types:" in out
    assert charts in out


def test_build_rag_retrieval_query_question_only():
    assert _build_rag_retrieval_query("list clusters", "") == "list clusters"


def test_build_rag_retrieval_query_charts_only():
    assert _build_rag_retrieval_query("", '"box"') == '"box"'
