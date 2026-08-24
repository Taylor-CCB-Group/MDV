# Server extensions

MDV database/catalog deployments select server extensions through the existing
`extensions` configuration array. Built-in extensions and separately installed
extensions share the same lifecycle and the same three-method interface.

## Enable an extension

Add the extension's unique ID to the deployment configuration:

```json
{
  "extensions": ["project_manager", "example_extension"]
}
```

Installation and enablement are separate. An installed external extension is not
imported unless its ID appears in this array. A configured extension that is
missing, duplicated, malformed, or fails during startup stops the application
with an `ExtensionError`.

The configured order is preserved after the always-on `ucsc_proxy` extension.
Each enabled extension is instantiated once per Flask application and that same
instance is reused for global and project behavior.

The host implementation keeps these responsibilities separate:

- `server_options.py` retains the built-in extension map, reads the existing
  configuration array, and passes the active instances to `MDVServerOptions`.
- `external_extensions.py` discovers and validates separately installed
  extensions. It is only asked to load names that the deployment enabled.
- `extension_navigation.py` validates optional navigation and exposes the
  catalog navigation endpoint.

## Publish an external extension

An external Python distribution advertises a zero-argument factory through the
`mdv.extensions` entry-point group:

```toml
[project.entry-points."mdv.extensions"]
example_extension = "example_package.mdv_extension:create_extension"
```

The returned object must conform to `MDVProjectServerExtension`:

```python
class ExampleExtension:
    def register_global_routes(self, app, config):
        ...

    def register_routes(self, project, project_bp):
        ...

    def mutate_state_json(self, state_json, project, app):
        ...


def create_extension():
    return ExampleExtension()
```

A global-only extension registers its Flask Blueprint in
`register_global_routes()` and can implement the two project methods as no-ops.
The external distribution may package and import its own dependencies normally.

Extension IDs must be unique. An external distribution cannot shadow an enabled
built-in extension, and multiple installed distributions cannot provide the same
enabled ID.

## Optional catalog navigation

An extension may contribute one catalog link without changing the required
three-method interface:

```python
from mdvtools.server_extension import ExtensionNavigation


class ExampleExtension:
    navigation = ExtensionNavigation(
        label="Example",
        url="/example/",
        requires_admin=True,
    )
```

The URL must be an application-relative path. When `requires_admin` is true, MDV
omits the link for anonymous and non-administrator users. Navigation visibility
does not protect the route itself.

Extension routes that require an administrator should use MDV's host-owned
authorization helper:

```python
from mdvtools.auth.authutils import admin_required


@blueprint.get("/api/example")
@admin_required
def example():
    ...
```

With authentication enabled, anonymous users are redirected to MDV login,
authenticated non-administrators receive HTTP 403, and administrators are
allowed. With authentication disabled, the helper allows local development
access.
