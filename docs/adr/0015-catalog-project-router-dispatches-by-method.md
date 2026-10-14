# The catalog's project router dispatches by method; routes without declared methods accept any

**Status:** accepted.

The catalog mounts every project's routes on `ProjectBlueprint`, a small regex router in
`project_router.py`, because Flask cannot add blueprints after the app has started serving. The
single-project server uses Flask's own routing through `SingleProjectShim`. `ProjectBlueprint`
stored one handler per path pattern and ignored the request method, so two routes on one path with
different methods, such as `POST /jobs` and `GET /jobs`, collapsed into one. In the catalog every
job submission ran the list handler and returned `200 []`.

We decided that `ProjectBlueprint` keeps every handler registered for a path and dispatches on the
request method: a route that declares `methods` handles only those methods, and when a path matches
but no handler accepts the method, the router returns 405, as Flask does. **A route registered
without `methods` keeps accepting any method.**

## Considered options

- **Match Flask exactly**, where a route without `methods` accepts only GET. The catalog and the
  single-project server would then route identically. We did not take it yet, because some routes
  are registered only in the catalog, for example in `project_manager_extension.py`, and have never
  run under Flask's rules. Any of them registered without `methods` and called with POST would start
  returning 405.
- **Merge `GET /jobs` and `POST /jobs` into one handler** that branches on `request.method`. This
  fixes the symptom for one route but leaves the router able to lose a handler silently the next
  time a path is shared.

## Consequences

- The two routers still differ for routes without `methods`: the catalog accepts any method and the
  single-project server accepts only GET. Matching Flask exactly is a follow-up that needs an audit
  of the routes that declare no methods.
- Within one method, routes are still tried in registration order, not by Werkzeug's specificity
  rules.
