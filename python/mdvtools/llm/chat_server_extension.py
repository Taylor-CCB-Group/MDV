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
    def register_routes(self, project: MDVProject, project_bp: ProjectBlueprintProtocol):
        # as well as registering routes, we might want to add something to be included in state.json?
        # or maybe we prefer to use websockets for this?
        # i.e. we always have a websocket connection to the client, 
        # and an extension can send a message to the client so it knows to update the UI.
        
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
        state_json["chat_enabled"] = chat_enabled

chat_extension = MDVProjectChatServerExtension()