"""
Dynamic Project Router for MDVTools

This module provides runtime project registration for Flask applications where
projects can be added/removed dynamically without restarting the server.

REFACTORING NOTES (for future work):
=====================================

ISSUES WITH CURRENT APPROACH:
- Custom regex-based routing reimplements Flask/Werkzeug functionality
- No type conversion for captured route parameters (int, float, etc.)
- Route matching order depends on dict iteration (predictable but not explicit)
- Global class-level state in ProjectBlueprint.blueprints
- Cannot use Flask's url_for() with blueprint-style endpoint names

REFACTORING OPTIONS TO CONSIDER:
1. Use werkzeug.routing.Map directly
   - Pro: Full Flask routing semantics, proper type conversion
   - Pro: Route priority/specificity rules built-in
   - Con: Need to rebuild Map on each project add/remove
   - Effort: Medium (2-3 days)

2. Flask PluggableView pattern
   - Pro: More Flask-native
   - Con: Still need custom dispatching for project-specific routes
   - Effort: Medium (2-3 days)

3. Dynamic Blueprint registration with app reload
   - Pro: True Flask Blueprints
   - Con: Requires app restart/reload on project changes
   - Con: May impact active connections/sessions
   - Effort: High (need graceful reload strategy)

4. Microservices approach
   - Pro: Each project could be independent service
   - Con: Major architectural change
   - Effort: Very High (weeks/months)

RECOMMENDED APPROACH:
If refactoring, start with Option 1 (werkzeug.routing.Map) as it provides
the best balance of proper Flask semantics without major architectural changes.

TESTS NEEDED BEFORE REFACTORING:
- Route pattern matching with all converter types
- before_request hook behavior
- Error handling and 404/500 responses
- Concurrent project add/remove operations
- Integration tests for auth hooks with route options
"""

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
IMPORTANT: This is a workaround for dynamic project registration at runtime.

WHY THIS EXISTS:
Flask Blueprints must be registered at application initialization time. However,
this application needs to add and remove MDVProjects dynamically while the server
is running. This class provides a Blueprint-like interface that supports runtime
route registration through a single catch-all route that dispatches to project-specific
handlers.

LIMITATIONS:
- Not a true Flask Blueprint - reimplements routing with regex pattern matching
- No support for url_for() with blueprint endpoints
- No blueprint-level error handlers or template folders
- Uses class-level dict for global state (ProjectBlueprint.blueprints)
- Route matching is custom and may not handle all Flask/Werkzeug edge cases

POTENTIAL REFACTORING OPTIONS (for future consideration):
1. Use Flask's PluggableView/MethodView for more Flask-native approach
2. Use werkzeug.routing.Map directly instead of custom regex matching
3. Consider if dynamic project registration is still required - if projects are
   mostly static, could move back to standard Blueprints with app reload on change
4. Investigate Flask extensions for dynamic routing (e.g., flask-classful)

CURRENT STATUS:
Works as a drop-in replacement for Blueprint in this codebase, but behavior differs
from real Blueprints. Use with caution and prefer standard Flask patterns where possible.
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
            r"""Convert Flask-style route rules to regex patterns.
            
            Supports Flask/Werkzeug converters:
            - <var> or <string:var>: matches text without slashes [^/]+
            - <path:var>: matches everything including slashes (.+)
            - <int:var>: matches integers (\d+)
            - <float:var>: matches floats (\d+\.?\d*)
            - <uuid:var>: matches UUIDs
            
            NOTE: This is a simplified implementation. For full Werkzeug routing semantics,
            consider using werkzeug.routing.Map directly in a future refactor.
            
            Current implementation details:
            - Query parameters are ignored (handled by trailing '/' in catch-all route)
            - Uses ^...$ anchors to match full path only
            - Captured groups are passed as positional args to view function
            """
            # Strategy: Replace <...> with placeholders, escape, then restore with regex patterns
            
            # Step 1: Extract all converter patterns and replace with unique placeholders
            converters = []
            placeholder_pattern = r'<([^>]+)>'
            
            def extract_converter(match):
                converters.append(match.group(1))
                return f'__CONVERTER_{len(converters) - 1}__'
            
            rule_with_placeholders = re.sub(placeholder_pattern, extract_converter, rule)
            
            # Step 2: Escape special regex characters (dots, etc.)
            rule_escaped = re.escape(rule_with_placeholders)
            
            # Step 3: Replace placeholders with actual regex patterns
            def get_converter_pattern(converter_spec):
                # Parse converter type and variable name
                if ':' in converter_spec:
                    converter, var_name = converter_spec.split(':', 1)
                else:
                    converter = 'string'  # Flask default
                    var_name = converter_spec
                
                # Map converters to regex patterns
                # These patterns align with Werkzeug's default converters
                patterns = {
                    'string': r'([^/]+)',    # Matches anything except / (Flask default)
                    'path': r'(.+)',         # Matches everything including /
                    'int': r'(\d+)',         # Matches integers only
                    'float': r'(\d+\.?\d*)', # Matches floats (simplified)
                    'uuid': r'([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})',
                }
                
                return patterns.get(converter, r'([^/]+)')  # Default to string if unknown
            
            # Replace each placeholder with its regex pattern
            rule_pattern = rule_escaped
            for i, converter_spec in enumerate(converters):
                placeholder = re.escape(f'__CONVERTER_{i}__')
                pattern = get_converter_pattern(converter_spec)
                rule_pattern = rule_pattern.replace(placeholder, pattern)
            
            rule_re = re.compile(f'^{rule_pattern}$')
            
            # Store route with handler and options (options may contain access_level, methods, etc.)
            self.routes[rule_re] = (f, options)
            return f

        return decorator

    def dispatch_request(self, subpath: str, project_id) -> Response:
        """Dispatch request to appropriate handler based on subpath pattern matching.
        
        Process:
        1. Normalize subpath to absolute path (strip query params)
        2. Try each registered route pattern in order
        3. On match, run before_request hooks (for auth/validation)
        4. If hooks pass, call view function with captured groups as args
        
        Example: '/tracks/<path:path>' compiles to regex '/tracks/(.+)'
                 Request '/tracks/mytrack/file.json' matches and calls handler('mytrack/file.json')
        
        NOTE: Routes are tried in iteration order (Python 3.7+ dict ordering).
        For more control over route priority, consider using werkzeug.routing.Map
        which supports priority/specificity rules.
        """
        # Parse URL to get clean path without query parameters
        subpath = f"/{urlparse(subpath).path}"
        
        for rule, (method, options) in self.routes.items():
            match = rule.match(subpath)
            if match:
                # Run all before_request functions (used for auth, logging, etc.)
                for func in self.before_request_funcs:
                    # Hooks receive project_id and route options (e.g., access_level)
                    rv = func(project_id=project_id, options=options)
                    if rv is not None:
                        # Hook returned a response - short-circuit and return it
                        # (typically used for auth failures, redirects, etc.)
                        return rv
                
                # All hooks passed - call the view function with captured groups
                # TODO: Consider type conversion here (int, float, etc.) based on converter
                # Currently just passes strings - handlers must convert themselves
                return method(*match.groups())
        
        # No route matched
        raise ValueError(f"no matching route for {subpath}")



class SingleProjectShim(ProjectBlueprintProtocol):
    """Adapter for single-project mode that uses standard Flask routing.
    
    When serving only one project (no dynamic multi-project support needed),
    this shim allows the same code to work by forwarding route() calls directly
    to Flask's native routing. No custom regex matching needed in this case.
    
    This is the "proper" Flask approach - use this pattern when possible.
    """
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
