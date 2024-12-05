from flask import Flask, Response, request, jsonify, session, redirect, url_for
from typing import Dict, Any, Callable, Tuple
import re
from urllib.parse import urlparse
from datetime import datetime, timedelta
import functools

"""
This should work as a drop-in replacement for `Blueprint` in the context
of a Flask app that can add (and remove) MDVProjects dynamically at runtime.
"""

class ProjectBlueprint:
    blueprints: Dict[str, "ProjectBlueprint"] = {}

    @staticmethod
    def register_app(app: Flask) -> None:
        @app.route(
            "/project/<project_id>/", defaults={"subpath": ""}, methods=["GET", "POST"]
        )
        @app.route("/project/<project_id>/<path:subpath>/", methods=["GET", "POST"])
        #incorporated below to resolve issue related to redirecting the request to http
        #@app.route("/project/<project_id>", defaults={'subpath': ''}, methods=["GET", "POST"])
        def project_route(project_id: str, subpath: str):
            """This generic route should call the appropriate method on the project with the given id.

            It will look up the project_id in ProjectBlueprint.blueprints and call the method with the given subpath.
            The ProjectBlueprint instance is responsible for routing the request to the correct method etc.
            """
            if project_id in ProjectBlueprint.blueprints:
                try:
                    return ProjectBlueprint.blueprints[project_id].dispatch_request(
                        subpath
                    )
                except Exception as e:
                    return {"status": "error", "message": str(e)}
            return {"status": "error", "message": "invalid project_id or method"}, 500

    def __init__(self, name: str, _ignored: str, url_prefix: str) -> None:
        self.name = name
        self.url_prefix = url_prefix  # i.e. /project/<project_id>/
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
            rule_re = re.compile(re.sub(r"<.*?>", "(.*)", f"^{rule}$"))
            self.routes[rule_re] = f
            return f

        return decorator

    def dispatch_request(self, subpath: str) -> Response:
        """We need to parse `subpath` so that we can route in a compatible way:
        Currently, we have regex patterns as keys in self.routes that match the subpath, with groups for
        any dynamic parts, so '/tracks/<path:path>' becomes '/tracks/(.*)'. If we get a request for '/tracks/mytrack',
        it will match the rule and call the method with 'mytrack' as the argument.
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
        subpath = f"/{urlparse(subpath).path}"
        for rule, method in self.routes.items():
            match = rule.match(subpath)
            if match:
                return method(*match.groups())
        raise ValueError(f"no matching route for {subpath}")
    
class ProjectBlueprint_v2:
    blueprints: Dict[str, "ProjectBlueprint_v2"] = {}
    # Class-level constant for the timestamp update interval
    TIMESTAMP_UPDATE_INTERVAL = timedelta(hours=1)

    @staticmethod
    def register_app(app: Flask) -> None:
        @app.route(
            "/project/<project_id>/", defaults={"subpath": ""}, methods=["GET", "POST"]
        )
        @app.route("/project/<project_id>/<path:subpath>/", methods=["GET", "POST"])
        #incorporated below to resolve issue related to redirecting the request to http
        #@app.route("/project/<project_id>", defaults={'subpath': ''}, methods=["GET", "POST"])
        def project_route(project_id: str, subpath: str):
            """This generic route should call the appropriate method on the project with the given id.

            It will look up the project_id in ProjectBlueprint_v2.blueprints and call the method with the given subpath.
            The ProjectBlueprint_v2 instance is responsible for routing the request to the correct method etc.
            """
            if not ProjectBlueprint_v2.is_authenticated():
                return redirect(url_for('login_dev'))
            
            if project_id in ProjectBlueprint_v2.blueprints:
                try:
                    return ProjectBlueprint_v2.blueprints[project_id].dispatch_request(
                        subpath,
                        project_id
                    )
                except Exception as e:
                    return {"status": "error", "message": str(e)}
            return {"status": "error", "message": "invalid project_id or method"}, 500
    
    @staticmethod
    def is_authenticated() -> bool:
        """Checks if the user is authenticated."""
        return 'token' in session

    def __init__(self, name: str, _ignored: str, url_prefix: str) -> None:
        self.name = name
        self.url_prefix = url_prefix  # i.e. /project/<project_id>/
        self.routes: Dict[re.Pattern[str], Tuple[Callable, Dict[str, Any]]] = {}
        ProjectBlueprint_v2.blueprints[name] = self

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
            rule_re = re.compile(re.sub(r"<.*?>", "(.*)", f"^{rule}$"))
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
        # find the method in self.routes that matches the subpath
        # e.g. '/tracks/mytrack' -> self.routes['/tracks/<path:path>']('mytrack')
        # /get_data -> self.routes['/get_data']()
        # '/datasources.json' -> self.routes['/<file>.json']('datasources')...
        # we need to parse the subpath to find the right method
        # as of now, we only have 0 or 1 <dynamic> parts in the rules
        # so we can match '/tracks/mytrack' to '/tracks/<path:path>'
        # and pass 'mytrack' to the method

        # find the item in self.routes that matches the subpath
        from mdvtools.dbutils.dbservice import ProjectService
        print("***********************************")
        print(subpath, project_id)
        subpath = f"/{urlparse(subpath).path}"
        for rule, (method, options) in self.routes.items():
            match = rule.match(subpath)
            if match:
                
                # Update the accessed timestamp only if the last update was more than TIMESTAMP_UPDATE_INTERVAL ago
                project = ProjectService.get_project_by_id(project_id)
                if project and (not project.accessed_timestamp or 
                    datetime.now() - project.accessed_timestamp > self.TIMESTAMP_UPDATE_INTERVAL):
                    print("****time interval greater than an hour ")
                    try:
                        ProjectService.set_project_accessed_timestamp(project_id)
                        print(f"In dispatch_request: Updated accessed timestamp for project ID {project_id}")
                    except Exception as e:
                        print(f"dispatch_request: Error updating accessed timestamp: {e}")
                        return jsonify({"status": "error", "message": "Failed to update project accessed timestamp"}), 500

                # first determine whether this is allowed
                # - rather than iterating over (rule, method), we might have
                # (rule, (method, permissionsFlags))
                # we can check the request.token (or whatever it is) here...
                print("options",options)
                # Check for access level only if specified in options
                if options and 'access_level' in options:
                    print("match", match)
                    #project_id = match.group(0).split('/')[2]  # Extract project_id from matched route
                    project = ProjectService.get_project_by_id(project_id)  # Fetch the project
                    
                    if project is None:
                        print("In dispatch_request: Error - project doesn't exist")
                        return jsonify({"status": "error", "message": f"Project with ID {project_id} not found"}), 404

                    required_access_level = options['access_level']  # Get required access level
                    print("required_access_level", required_access_level)
                    if required_access_level == 'editable':
                        print("required_access_level is editable, fetched project is ", project)
                        if project.access_level != 'editable':
                            return jsonify({"status": "error", "message": "This project is read-only and cannot be modified."}), 403
                
                return method(*match.groups())
            
        raise ValueError(f"no matching route for {subpath}")



class SingleProjectShim:
    def __init__(self, app: Flask) -> None:
        self.app = app

    def route(self, rule: str, **options: Any) -> Callable:
        access_level = options.pop("access_level", None)  # Remove access_level if present

        def decorator(func: Callable) -> Callable:
            # Define a wrapper around the original view function
            @functools.wraps(func)  # Preserve the original function’s name and docstring
            def wrapped_func(*args, **kwargs):
                # Process access_level or other logic here if needed
                if access_level:
                    print(f"Access level required: {access_level}")
                # Call the original function
                return func(*args, **kwargs)

            # Set a unique endpoint name for each wrapped function using func's original name
            endpoint = options.pop("endpoint", func.__name__)
            # Register the route with Flask, using the unique endpoint name
            return self.app.route(rule, endpoint=endpoint, **options)(wrapped_func)

        return decorator
