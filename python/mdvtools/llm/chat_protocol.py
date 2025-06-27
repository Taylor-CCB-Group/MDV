from typing import Callable, Optional, Protocol, Any, TypedDict
from mdvtools.mdvproject import MDVProject


class AskQuestionResult(TypedDict):
    code: Optional[str]
    view_name: Optional[str]
    error: bool
    message: str


class ProjectChatProtocol(Protocol):
    def __init__(self, project: MDVProject): ...
    def log(self, msg: str, *args: Any) -> None: ...

    def ask_question(
            self, 
            question: str, 
            id: str, 
            conversation_id: str, 
            room: str,
            handle_error: Callable[[str], None],
            ) -> AskQuestionResult: ...
    welcome: str
    init_error: Optional[bool]
    error_message: Optional[str]

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

        def log(self, msg, *args):
            print("[dummy bot]: " + (msg % args if args else msg))

        def ask_question(self, question, id, conversation_id, room, handle_error):
            handle_error(f"Sorry, I can't help you right now\n\n{msg}")
            return AskQuestionResult(
                code=None,
                view_name=None,
                error=True,
                message=f"Sorry, I can't help you right now\n\n{msg}"
            )

# Tell the type checker we’re exposing a ProjectChat conforming to the protocol
ProjectChat: type[ProjectChatProtocol] = _ProjectChat
print(f">>> chat_enabled: {chat_enabled}")
