from flask import Flask
from flask.typing import ResponseReturnValue as Response
from typing import Dict, Any, Callable, Tuple, List
import re
from urllib.parse import urlparse
import functools
from typing import Protocol
from mdvtools.logging_config import get_logger

logger = get_logger(__name__)

"""
This should work as a drop-in replacement for `Blueprint` in the context
of a Flask app that can add (and remove) MDVProjects dynamically at runtime.
"""
class ProjectBlueprintProtocol(Protocol):
    def route(self, rule: str, **options: Any) -> Callable:
        ...
    
    def before_request(self, f: Callable) -> Callable:
        ...

class ProjectBlueprint(ProjectBlueprintProtocol):
    blueprints: Dict[str, "ProjectBlueprint"] = {}

    @staticmethod
    def register_app(app: Flask) -> None:
        @app.route(
            "/project/<project_id>/", defaults={"subpath": ""}, methods=["GET", "POST", "PATCH"]
        )
        @app.route("/project/<project_id>/<path:subpath>/", methods=["GET", "POST", "PATCH"])
        def project_route(project_id: str, subpath: str):
            """This generic route should call the appropriate method on the project with the given id.

            It will look up the project_id in ProjectBlueprint.blueprints and call the method with the given subpath.
            The ProjectBlueprint instance is responsible for routing the request to the correct method etc.
            """
            if project_id in ProjectBlueprint.blueprints:
                try:
                    return ProjectBlueprint.blueprints[project_id].dispatch_request(
                        subpath, project_id
                    )
                except Exception as e:
                    logger.exception(f"Error dispatching request for project {project_id} on subpath {subpath}")
                    return {"status": "error", "message": str(e)}, 500
            return {"status": "error", "message": "invalid project_id or method"}, 500

    def __init__(self, name: str, _ignored: str, url_prefix: str) -> None:
        self.name = name
        self.url_prefix = url_prefix  # i.e. /project/<project_id>/
        self.routes: Dict[re.Pattern[str], Tuple[Callable, Dict[str, Any]]] = {}
        self.before_request_funcs: List[Callable] = []
        ProjectBlueprint.blueprints[name] = self

    def before_request(self, f: Callable) -> Callable:
        """Register a function to run before each request for this blueprint."""
        self.before_request_funcs.append(f)
        return f

    def route(self, rule: str, **options: Any) -> Callable:
        def decorator(f: Callable) -> Callable:
            """As of this writing, the rules in server.py all have one or zero dynamic parts.
            i.e. all of our view_func have 0 or 1 arguments.
            We use route keys which are regex patterns to match the subpath and parse out the dynamic parts.
            We should also ignore any query parameters in the subpath (trailing '/' in generic overall route does this).
            We have a further issue with '/' matching all routes, so we use f'^{rule}$' to match the full path.
            """
            # would be nice to re-use logic from Flask's routing, but it's a bit complex to use out of context
            # from werkzeug.routing import Rule
            # Escape literal '.' in the rule to '\.' so that dots match dots, not any character
            rule_escaped = rule.replace('.', r'\.')
            rule_re = re.compile(re.sub(r"<.*?>", "(.*)", f"^{rule_escaped}$"))
            # rather than just add the callable f, we add (f, permissionsFlags)
            # where permissionsFlags are something passed in options
            self.routes[rule_re] = (f, options)
            return f

        return decorator

    def dispatch_request(self, subpath: str, project_id) -> Response:
        """We need to parse `subpath` so that we can route in a compatible way:
        Currently, we have regex patterns as keys in self.routes that match the subpath, with groups for
        any dynamic parts, so '/tracks/<path:path>' becomes '/tracks/(.*)'. If we get a request for '/tracks/mytrack',
        it will match the rule and call the method with 'mytrack' as the argument.
        """
        subpath = f"/{urlparse(subpath).path}"
        for rule, (method, options) in self.routes.items():
            match = rule.match(subpath)
            if match:
                # Run all before_request functions. These are responsible for auth, etc.
                for func in self.before_request_funcs:
                    # The hook needs access to project_id and the route options
                    rv = func(project_id=project_id, options=options)
                    if rv is not None:
                        # If a hook returns a value, it becomes the response, short-circuiting the request
                        return rv
                
                # If all hooks passed (returned None), call the original view function
                return method(*match.groups())
            
        raise ValueError(f"no matching route for {subpath}")



class SingleProjectShim(ProjectBlueprintProtocol):
    def __init__(self, app: Flask) -> None:
        self.app = app

    def before_request(self, f: Callable) -> Callable:
        # This is a no-op for the single project shim.
        # It's here to satisfy the ProjectBlueprintProtocol.
        # logger.warning("before_request is not implemented for SingleProjectShim and has no effect.")
        return f

    def route(self, rule: str, **options: Any) -> Callable:
        access_level = options.pop("access_level", None)  # Remove access_level if present

        def decorator(func: Callable) -> Callable:
            # Define a wrapper around the original view function
            @functools.wraps(func)  # Preserve the original function's name and docstring
            def wrapped_func(*args, **kwargs):
                # Process access_level or other logic here if needed
                if access_level:
                    logger.info(f"Access level required: {access_level}")
                # Call the original function
                return func(*args, **kwargs)

            # Set a unique endpoint name for each wrapped function using func's original name
            endpoint = options.pop("endpoint", func.__name__)
            # Register the route with Flask, using the unique endpoint name
            return self.app.route(rule, endpoint=endpoint, **options)(wrapped_func)

        return decorator
