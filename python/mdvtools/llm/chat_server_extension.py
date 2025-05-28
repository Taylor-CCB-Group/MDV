from typing import Optional, Protocol
from mdvtools.llm.chat_protocol import (
  ProjectChat,
  ProjectChatProtocol, 
  chat_enabled
)
from mdvtools.llm.code_manipulation import parse_view_name
from mdvtools.mdvproject import MDVProject
from flask import request
from mdvtools.project_router import ProjectBlueprintProtocol
from mdvtools.websocket import socketio
# from mdvtools.dbutils.config import config
import os

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
    def register_routes(self, project: MDVProject, blueprint: ProjectBlueprintProtocol):
        """
        Assign any extra `/project/<project_id>/<path>` routes to the blueprint for this project instance.
        """
        ...
    def mutate_state_json(self, state_json: dict, project: MDVProject):
        """
        Mutate the state.json before returning it as a request response,
        e.g. to add information about the extension.
        """
        ...


class MDVProjectChatServerExtension(MDVProjectServerExtension):
    # as well as registering routes and mutating the state.json,
    # we might describe websocket routes here, and how we control auth with that
    def register_routes(self, project: MDVProject, project_bp: ProjectBlueprintProtocol):
        # @socketio.on("connect", namespace=f"/project/{project.id}")
        # def chat_connect():
        #     """
        #     Handle WebSocket connections for the chat extension.
        #     This could be used to initialize a chat session or perform other setup tasks.
        #     We should check authentication here,
        #     """
            

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
            try:
                if bot is None:
                    bot = ProjectChat(project)
                # we need to know view_name as well as message - but also maybe there won't be one, if there's an error etc.
                # probably want to change the return type of this function, but for now we do some string parsing here.
                final_code = bot.ask_question(message, id)
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
            return {"message": f"bleep bloop I'm a robot, you said: {message}"}
    def mutate_state_json(self, state_json: dict, project: MDVProject):
        """
        Mutate the state.json before returning it as a request response,
        in this case to add information about the chat extension.
        This should also reflect what the user is authorized to do,
        e.g. whether they can use the chat functionality or not.
        Currently, all chat related routes will require the user to have 'editable' access level,
        although we don't currently signal to the front-end whether the user has this access level or not.
        """
        was_enabled = state_json.get("chat_enabled", False)
        state_json["chat_enabled"] = chat_enabled and was_enabled
        # I think we want to maybe get something from a global config for e.g. server_root="https://server.example/container_path"
        # for now, we can probably patch this in the front-end code without passing this back.
        # state_json["chat_route"] = f"/project/{project.id}/chat"
        # hard-coded default for demo purposes...
        state_json["socketio_route"] = os.getenv("SOCKETIO_ROUTE", "https://bia.cmd.ox.ac.uk/socket.io")

chat_extension = MDVProjectChatServerExtension()