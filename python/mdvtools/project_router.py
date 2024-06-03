from flask import Flask, Response
from typing import Dict, Any, Callable
import re
from urllib.parse import urlparse

"""
This should work as a drop-in replacement for `Blueprint` in the context
of a Flask app that can add (and remove) MDVProjects dynamically at runtime.
"""
class ProjectBlueprint:
    blueprints: Dict[str, "ProjectBlueprint"] = {}
    @staticmethod
    def register_app(app: Flask) -> None:
        # return
        print("registering app")
        @app.route("/project/<project_id>/", defaults={'subpath': ''}, methods=["GET", "POST"])
        @app.route("/project/<project_id>/<subpath>/", methods=["GET", "POST"])
        def project_route(project_id: str, subpath: str):
            """This generic route should call the appropriate method on the project with the given id.
            
            It will look up the project_id in ProjectBlueprint.blueprints and call the method with the given subpath.
            The ProjectBlueprint instance is responsible for routing the request to the correct method etc.
            """
            print(f'project_route: {project_id}, {subpath}')
            if project_id in ProjectBlueprint.blueprints:
                try:
                    print(f'found project {project_id}, dispatching request to {subpath}...')
                    return ProjectBlueprint.blueprints[project_id].dispatch_request(subpath)
                except Exception as e:
                    return {"status": "error", "message": str(e)}
            return {"status": "error", "message": "invalid project_id or method"}, 500

    def __init__(self, name: str, _ignored: str, url_prefix: str) -> None:
        self.name = name
        self.url_prefix = url_prefix # i.e. /project/<project_id>/
        self.routes: Dict[re.Pattern[str], Callable] = {}
        ProjectBlueprint.blueprints[name] = self

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
            rule_re = re.compile(re.sub(r'<.*?>', '(.*)', f'^{rule}$'))
            self.routes[rule_re] = f
            return f

        return decorator
    
    def dispatch_request(self, subpath: str) -> Response:
        """We need to parse `subpath` so that we can route in a compatible way:
        If we have a route like `/tracks/<path:path>`, we'll receive `path` as a single string.
        We need to store the method in blueprint.routes so that we can find it when user requests
        `/project/<project_id>/tracks/mytrack`... right now, we'll try to look for '/tracks/mytrack'
        We need to split it into a list of strings to pass to the method.
        """
        # find the method in self.routes that matches the subpath
        # e.g. '/tracks/mytrack' -> self.routes['/tracks/<path:path>']('mytrack')
        # /get_data -> self.routes['/get_data']()
        # '/datasources.json' -> self.routes['/<file>.json']('datasources')...
        # we need to parse the subpath to find the right method
        # as of now, we only have 0 or 1 <dynamic> parts in the rules
        # so we can match '/tracks/mytrack' to '/tracks/<path:path>'
        # and pass 'mytrack' to the method
        
        # find the item in self.routes that matches the subpath
        subpath = f'/{urlparse(subpath).path}'
        for rule, method in self.routes.items():
            match = rule.match(subpath)
            if match:
                print(f'matching route {rule.pattern} to {subpath}...')
                return method(*match.groups())
        raise ValueError(f"no matching route for {subpath}")