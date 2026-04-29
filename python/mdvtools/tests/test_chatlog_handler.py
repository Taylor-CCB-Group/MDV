import logging

from mdvtools.llm.chatlog import LangchainLoggingHandler


class _DummyResponse:
    def __init__(self) -> None:
        self.generations = [[]]


def test_langchain_logging_handler_avoids_info_noise(caplog):
    logger = logging.getLogger("test_langchain_logging_handler")
    logger.handlers = []
    logger.setLevel(logging.INFO)
    handler = LangchainLoggingHandler(logger)

    with caplog.at_level(logging.INFO):
        handler.on_chat_model_start({}, [])
        handler.on_llm_end(_DummyResponse())
        handler.on_chain_start({"name": "demo"}, {})
        handler.on_chain_end({})

    assert "Chat model started" not in caplog.text
    assert "LLMResult" not in caplog.text
    assert "Chain ended" not in caplog.text
