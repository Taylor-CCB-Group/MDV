import logging

from langchain_core.outputs import LLMResult

from mdvtools.llm.chatlog import LangchainLoggingHandler


def test_langchain_logging_handler_avoids_info_noise(caplog):
    logger = logging.getLogger("test_langchain_logging_handler")
    logger.handlers = []
    logger.setLevel(logging.INFO)
    handler = LangchainLoggingHandler(logger)

    with caplog.at_level(logging.INFO):
        handler.on_chat_model_start({}, [])
        handler.on_llm_end(LLMResult(generations=[[]]))
        handler.on_chain_start({"name": "demo"}, {})
        handler.on_chain_end({})

    assert "Chat model started" not in caplog.text
    assert "LLMResult" not in caplog.text
    assert "Chain ended" not in caplog.text
