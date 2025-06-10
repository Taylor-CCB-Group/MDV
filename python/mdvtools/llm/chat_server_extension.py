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
from flask import Flask, request, session, current_app

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
            

        bot: Optional[ProjectChatProtocol] = None
        @project_bp.route("/chat_init", access_level='editable', methods=["POST"])
        def chat_init():
            nonlocal bot
            if bot is None:
                bot = ProjectChat(project)
            return {"message": bot.welcome}
        @project_bp.route("/chat", access_level='editable', methods=["POST"])
        def chat():
            nonlocal bot
            if not request.json:
                return {"error": "No JSON data in request"}, 500
            message = request.json.get("message")
            id = request.json.get("id")
            if not message or not id:
                return {"error": "Missing 'message' or 'id' in request JSON"}, 400
            # todo - consider having a socket.io event for this, rather than a REST endpoint.
            # this would mean that we could use request.sid to track the chat session
            conversation_id = request.json.get("conversation_id")
            try:
                if bot is None:
                    bot = ProjectChat(project)
                # we need to know view_name as well as message - but also maybe there won't be one, if there's an error etc.
                # probably want to change the return type of this function, but for now we do some string parsing here.
                final_code = bot.ask_question(message, id, conversation_id)
                try:
                    # this can give a confusing error if we don't explicitly catch it...
                    view_name = parse_view_name(final_code)
                    if view_name is None:
                        raise Exception(final_code)
                    bot.log(f"view_name: {view_name}")
                    return {"message": final_code, "view": view_name, "id": id}
                except Exception as e:
                    bot.log(f"final_code returned by bot.ask_question is bad, probably an earlier error: {e}")
                    # final_code is probably an error message, at this point.
                    return {"message": final_code}
            except Exception as e:
                print(e)
                return {"message": str(e)}
    
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