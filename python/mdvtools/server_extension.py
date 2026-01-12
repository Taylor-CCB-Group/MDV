from __future__ import annotations
from typing import Protocol, TYPE_CHECKING, List, Optional
from flask import Flask
from dataclasses import dataclass, field

if TYPE_CHECKING:
    from mdvtools.project_router import ProjectBlueprintProtocol
    from mdvtools.mdvproject import MDVProject


class MDVProjectServerExtension(Protocol):
    """
    A protocol for server extensions that can be used with MDV projects.

    We might use this for blocks of other functionality that aren't totally core mdv functionality,
    integrating other services/libraries/functionality.
    We might also re-arrange so that some things like the add_anndata routes are moved into an extension.

    Maybe rather  than pass a Flask app to `server.py`, we pass something representing MDV app configuration,
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
    def register_global_routes(self, app: Flask, config: dict):
        """
        Register global routes that apply to the entire application.
        """
        ...


@dataclass
class MDVServerOptions:
    """
    Configuration options for the MDV server.
    """
    open_browser: bool = True
    port: int = 5050
    websocket: bool = True
    use_reloader: bool = False
    app: Optional[Flask] = None
    backend_db: bool = False
    extensions: List[MDVProjectServerExtension] = field(default_factory=list)

