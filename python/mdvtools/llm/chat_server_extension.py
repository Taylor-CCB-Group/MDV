from typing import Optional, Protocol
from mdvtools.llm.chat_protocol import (
  ChatRequest,
  ProjectChat,
  ProjectChatProtocol, 
  chat_enabled
)
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprintProtocol
# from mdvtools.dbutils.config import config
from flask import Flask, request
from flask_socketio import join_room, leave_room
from mdvtools.logging_config import get_logger
from mdvtools.llm.chatlog import log_chat_item

logger = get_logger(__name__)
logger.info("chat server extension loading...")


class MDVProjectServerExtension(Protocol):
    """
    A protocol for server extensions that can be used with MDV projects.

    We might use this for blocks of other functionality that aren't totally core mdv functionality,
    integrating other services/libraries/functionality.
    We might also re-arrange so that some things like the add_anndata routes are moved into an extension.

    Maybe rather than pass a Flask app to `server.py`, we pass something representing MDV app configuration,
    including a Flask app, these extensions, auth provider etc...
    Flask becomes an implementation detail that we abstract away somewhat.
    """
    def register_routes(self, project: MDVProject, project_bp: ProjectBlueprintProtocol):
        """
        Assign any extra `/project/<project_id>/<path>` routes to the blueprint for this project instance.
        """
        ...
    def mutate_state_json(self, state_json: dict, project: MDVProject, app: Flask):
        """
        Mutate the state.json before returning it as a request response,
        e.g. to add information about the extension.

        Don't really want to pass flask app here, doing so for now to allow access to config.
        """
        ...


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
                try:
                    bot = ProjectChat(project)
                except Exception as e:
                    print(f"ERROR: {str(e)[:500]}")
                    return {"message": f"ERROR: {str(e)[:500]}", "error": True}
            if bot.init_error:
                # Log and return the error
                error_msg = bot.error_message or "An unknown error occurred"
                bot.log(f"ERROR: {error_msg}")
                return {"message": f"ERROR: {error_msg}", "error": True}
            return {"message": bot.welcome}

        @socketio.on("chat_request", namespace=f"/project/{project.id}")
        def chat(data):
            nonlocal bot
            sid = request.sid # type: ignore - flask_socketio.request.sid
            # todo - auth (via custom decorator? see comments above)

            message = data.get("message")
            id = data.get("id")
            room = f"{sid}_{id}"
            join_room(room)
            def handle_error(error: str):
                #! todo - frontend should handle this, and show an error message to the user.
                # also, this method may be augmented so that it also logs to the chat log,
                # and pass this handler to the bot.ask_question method.
                # Log error to chat_log.json
                log_chat_item(project, message or '', None, '', error, conversation_id, None, error=True)
                socketio.emit(
                    "chat_error",
                    {"message": error},
                    namespace=f"/project/{project.id}",
                    to=room
                )
            if not message or not id:
                handle_error("Missing 'message' or 'id' in request JSON")
                leave_room(room)
                return
            conversation_id = data.get("conversation_id")
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
                # print(f"ERROR: {str(e)[:500]}")
                handle_error(f"ERROR: {str(e)[:500]}")
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
        

chat_extension = MDVProjectChatServerExtension()