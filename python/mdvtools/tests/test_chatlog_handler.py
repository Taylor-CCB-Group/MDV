import logging

from langchain_core.outputs import LLMResult

from mdvtools.llm.chatlog import LangchainLoggingHandler


def test_chat_socket_api_log_file_only_skips_socket_handler():
    from unittest.mock import MagicMock, patch

    from mdvtools.llm.chatlog import ChatSocketAPI, ChatSocketIOHandler

    emitted: list[str] = []

    class FakeSocketIO:
        def emit(self, event, msg, namespace=None, to=None):
            emitted.append(msg)

    project = MagicMock()
    project.id = "proj1"
    project.dir = "/tmp/chatlog_test_project"

    with patch("mdvtools.websocket.socketio", FakeSocketIO()):
        api = ChatSocketAPI(project, "req1", "room1", "conv1")
        socket_handlers = [
            h for h in api.logger.handlers if isinstance(h, ChatSocketIOHandler)
        ]
        assert socket_handlers, "expected a ChatSocketIOHandler on the request logger"

        api.log_file_only("model_id=openai:chat:gpt-4.1 provider=openai")
        api.logger.info("visible on socket")

        assert emitted == ["visible on socket"]


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
