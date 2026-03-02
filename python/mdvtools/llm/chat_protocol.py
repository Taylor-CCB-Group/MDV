from typing import Optional, Protocol, Any, TypedDict, Union
from mdvtools.mdvproject import MDVProject


class AskQuestionResult(TypedDict):
    code: Optional[str]
    view_name: Optional[str]
    error: bool
    message: str


class HandleError(Protocol):
    def __call__(self, error: Union[str, Exception], *, extra_metadata: Optional[dict] = None) -> None:
        ...


class ChatRequest(TypedDict):
    message: str
    id: str
    conversation_id: str
    room: str
    handle_error: HandleError

class ProjectChatProtocol(Protocol):
    def __init__(self, project: MDVProject): ...
    def log(self, msg: str, *args: Any) -> None: ...
    def get_suggested_questions(self) -> list[str]: ...
    def ask_question(
            self, 
            chat_request: ChatRequest,
            ) -> AskQuestionResult: ...
    welcome: str
    init_error: Optional[bool]
    error_message: Optional[str]

chat_enabled = True # pending server extension config mechanism
try:
    from mdvtools.llm.langchain_mdv import ProjectChat as _ProjectChat # type: ignore
    # raise Exception("test error in import langchain_mdv")
except Exception as e:
    msg = str(e)
    class _ProjectChat(ProjectChatProtocol):
        def __init__(self, project: MDVProject):
            self.project = project
            self.init_error = True
            self.welcome = (
                "There was a problem initializing the chatbot:\n\n"
                f"{msg}"
            )
            self.error_message = msg

        def log(self, msg, *args):
            print("[dummy bot]: " + (msg % args if args else msg))

        def get_suggested_questions(self):
            return []

        def ask_question(self, chat_request: ChatRequest):
            handle_error = chat_request["handle_error"]
            msg = self.welcome
            handle_error(f"## Sorry, I can't help you right now...\n\n{msg}")
            return AskQuestionResult(
                code=None,
                view_name=None,
                error=True,
                message=f"Sorry, I can't help you right now\n\n{msg}"
            )

# Tell the type checker we’re exposing a ProjectChat conforming to the protocol
ProjectChat: type[ProjectChatProtocol] = _ProjectChat
