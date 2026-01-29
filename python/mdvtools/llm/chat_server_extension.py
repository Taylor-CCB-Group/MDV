import traceback
from typing import Optional, Union
from mdvtools.benchmarking.memory import log_memory_usage
from mdvtools.llm.chat_protocol import (
  ChatRequest,
  ProjectChat,
  ProjectChatProtocol, 
  chat_enabled
)
from mdvtools.markdown_utils import create_project_markdown, create_error_markdown
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprintProtocol
# from mdvtools.dbutils.config import config
from flask import Flask, request
from flask_socketio import join_room, leave_room
from mdvtools.logging_config import get_logger
from mdvtools.llm.chatlog import log_chat_item
from mdvtools.server_extension import MDVProjectServerExtension
from mdvtools.auth.authutils import is_authenticated

logger = get_logger(__name__)
logger.info("chat server extension loading...")



class MDVProjectChatServerExtension(MDVProjectServerExtension):
    # as well as registering routes and mutating the state.json,
    # we might describe websocket routes here, and how we control auth with that

    def register_routes(self, project: MDVProject, project_bp: ProjectBlueprintProtocol):
        from mdvtools.websocket import socketio
        if socketio is None:
            raise Exception("socketio is not initialized")
        @socketio.on("connect", namespace=f"/project/{project.id}")
        def chat_connect():
            """
            Handle WebSocket connections for the chat extension.
            This could be used to initialize a chat session or perform other setup tasks.
            We should check authentication here,
            **nb the coupling of sockets and "chat" should be removed.**
            """
            if not is_authenticated():
                return False
            # todo - check access level, whether chat is enabled etc
            # at the time you connect, flask_socketio session by default is "forked from http session"
            # if we made SocketIO instance with manage_session=False, then they would reference the same session.
            # Not sure how this relates to how we handle auth, access_level etc.
            # https://blog.miguelgrinberg.com/post/flask-socketio-and-the-user-session
            # https://flask-socketio.readthedocs.io/en/latest/implementation_notes.html#using-flask-login-with-flask-socketio
            

        bot: Optional[ProjectChatProtocol] = None
        @project_bp.route("/chat_init", access_level='editable', methods=["POST"])
        def chat_init():
            nonlocal bot
            if bot is None:
                log_memory_usage("before bot init")
                try:
                    bot = ProjectChat(project)
                    log_memory_usage("after bot init")
                except Exception as e:
                    error_message = f"{str(e)[:500]}\n\n{traceback.format_exc()}"
                    print(f"ERROR: {error_message}")
                    return {"message": f"ERROR: {error_message}", "error": True}
            if bot.init_error:
                # Log and return the error
                error_msg = bot.error_message or "An unknown error occurred"
                bot.log(f"ERROR: {error_msg}")
                return {"message": f"ERROR: {error_msg}", "error": True}
            
            markdown = create_project_markdown(project)
                
            bot.log(f"Markdown: {markdown}")
            detailed_message = bot.welcome + markdown
            # we want some suggested questions here.
            suggested_questions = bot.get_suggested_questions()
            logger.info(f"Suggested questions: {suggested_questions}")
            return {"message": detailed_message, "suggested_questions": suggested_questions}

        @socketio.on("chat_request", namespace=f"/project/{project.id}")
        def chat(data):
            if not is_authenticated():
                return

            nonlocal bot
            sid = request.sid # type: ignore - flask_socketio.request.sid
            # todo - auth (via custom decorator? see comments above)

            message = data.get("message")
            id = data.get("id")
            conversation_id = data.get("conversation_id")
            room = f"{sid}_{id}"
            join_room(room)
            def handle_error(error: Union[str, Exception], *, extra_metadata: Optional[dict] = None):
                if isinstance(error, Exception):
                    error_message = str(error)
                    traceback_str = traceback.format_exc()
                else:
                    error_message = error
                    traceback_str = None

                markdown = create_error_markdown(error_message, traceback_str, extra_metadata)
                logger.error(f"Chat error: {markdown}")
                log_chat_item(project, message or '', None, '', markdown, conversation_id, None, None, error=True)
                socketio.emit(
                    "chat_error",
                    {"message": markdown},
                    namespace=f"/project/{project.id}",
                    to=room
                )
            if not message or not id:
                handle_error("Missing 'message' or 'id' in request JSON")
                leave_room(room)
                return
            chat_request = ChatRequest(
                message=message,
                id=id,
                conversation_id=conversation_id,
                room=room,
                handle_error=handle_error
            )
            try:
                if bot is None:
                    # todo - allow this to be freed at some point if we're not using it anymore.
                    bot = ProjectChat(project)
                # we need to know view_name as well as message - but also maybe there won't be one, if there's an error etc.
                # probably want to change the return type of this function, but for now we do some string parsing here.
                result = bot.ask_question(chat_request)
                log_memory_usage("after bot ask_question")
                if result["error"]:
                    socketio.emit(
                        "chat_error", 
                        {"message": result["message"], "error": True, "id": id}, 
                        namespace=f"/project/{project.id}", 
                        to=room
                    )
                    leave_room(room)
                    return
                else:
                    socketio.emit(
                        "chat_response", 
                        {"message": result["code"], "view": result["view_name"], "id": id}, 
                        namespace=f"/project/{project.id}", 
                        to=room
                    )
                    leave_room(room)
                    return

            except Exception as e:
                # Log error to chat_log.json
                error_message = str(e)
                tb_str = traceback.format_exc()
                print(f"ERROR: {error_message}\n\n{tb_str}")
                handle_error(e)
                leave_room(room)
                return
            
    def mutate_state_json(self, state_json: dict, project: MDVProject, app: Flask):
        """
        Mutate the state.json before returning it as a request response,
        in this case to add information about the chat extension.
        This should also reflect what the user is authorized to do,
        e.g. whether they can use the chat functionality or not.
        Currently, all chat related routes will require the user to have 'editable' access level,
        although we don't currently signal to the front-end whether the user has this access level or not.
        """
        auth_enabled = app.config.get('ENABLE_AUTH', False)
        # allow chat functionality by default, unless auth is enabled
        was_enabled = state_json.get("chat_enabled", not auth_enabled)
        state_json["chat_enabled"] = chat_enabled and was_enabled
        
    def get_session_config(self):
        return {
            "chat_enabled": chat_enabled
        }
