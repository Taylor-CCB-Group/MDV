from typing import Optional, Protocol, cast
from mdvtools.llm.chat_protocol import (
  ProjectChat,
  ProjectChatProtocol, 
  chat_enabled
)
from mdvtools.llm.code_manipulation import parse_view_name
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprintProtocol
# from mdvtools.dbutils.config import config
from flask import Flask, request
from flask_socketio import join_room, leave_room

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
                bot = ProjectChat(project)
            return {"message": bot.welcome}

        @socketio.on("chat_request", namespace=f"/project/{project.id}")
        def chat(data):
            nonlocal bot
            sid = request.sid # type: ignore
            # todo - auth

            message = data.get("message")
            id = data.get("id")
            room = f"{sid}_{id}"
            join_room(room)
            def send_error(error: str):
                #! todo - frontend should handle this, and show an error message to the user.
                # also, this method may be augmented so that it also logs to the chat log,
                # and pass this handler to the bot.ask_question method.
                socketio.emit(
                    "chat_error",
                    {"message": error},
                    namespace=f"/project/{project.id}",
                    to=room
                )
            if not message or not id:
                send_error("Missing 'message' or 'id' in request JSON")
                leave_room(room)
                return
            conversation_id = data.get("conversation_id")
            try:
                if bot is None:
                    # todo - allow this to be freed at some point if we're not using it anymore.
                    bot = ProjectChat(project)
                # we need to know view_name as well as message - but also maybe there won't be one, if there's an error etc.
                # probably want to change the return type of this function, but for now we do some string parsing here.
                final_code = bot.ask_question(message, id, conversation_id, room)
                try:
                    # this can give a confusing error if we don't explicitly catch it...
                    view_name = parse_view_name(final_code)
                    if view_name is None:
                        raise Exception(final_code)
                    # oops - this log is no longer the right thing to do...
                    bot.log(f"view_name: {view_name}")
                    socketio.emit("chat_response", {"message": final_code, "view": view_name, "id": id}, namespace=f"/project/{project.id}", to=sid)
                    # return {"message": final_code, "view": view_name, "id": id}
                except Exception as e:
                    # oops - this log is no longer the right thing to do...
                    bot.log(f"final_code returned by bot.ask_question is bad, probably an earlier error: {e}")
                    send_error(str(e))
                    leave_room(room)
                    # final_code is probably an error message, at this point.
                    return
            except Exception as e:
                print(e)
                send_error(str(e))
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