from typing import Protocol, Any
from mdvtools.mdvproject import MDVProject

class ProjectChatProtocol(Protocol):
    def __init__(self, project: MDVProject): ...
    def log(self, msg: str, *args: Any, **kwargs: Any) -> None: ...
    def ask_question(self, message: str, id: str) -> str: ...
    welcome: str

chat_enabled = False
try:
    from mdvtools.llm.langchain_mdv import ProjectChat as _ProjectChat # type: ignore
    chat_enabled = True
except Exception as e:
    msg = str(e)
    class _ProjectChat(ProjectChatProtocol):
        def __init__(self, project: MDVProject):
            self.project = project
            self.welcome = (
                "There was a problem initializing the chatbot\n\n"
                f"{msg}"
            )

        def log(self, msg, *args, **kwargs):
            print("[dummy bot]: " + (msg % args if args else msg))

        def ask_question(self, message, id):
            return (
                "Sorry, I can't help you right now\n\n"
                f"{msg}"
            )

# Tell the type checker we’re exposing a ProjectChat conforming to the protocol
ProjectChat: type[ProjectChatProtocol] = _ProjectChat
print(f">>> chat_enabled: {chat_enabled}")
