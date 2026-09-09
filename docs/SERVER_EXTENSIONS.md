# Server extensions

MDV database/catalog deployments select server extensions through the existing
`extensions` configuration array. Built-in extensions and separately installed
extensions share the same lifecycle and the same three-method interface.

The design rationale, alternatives and consequences are recorded in
[ADR-0003](adr/0003-discover-installed-server-extensions.md). This document is
the operational reference for configuring and implementing an extension.

## Design summary

MDV separates three concerns that are easy to conflate:

- **Packaging:** a Python distribution contains extension code and dependencies.
- **Discovery:** an installed distribution advertises an extension factory through
  the `mdv.extensions` entry-point group.
- **Enablement:** a deployment explicitly lists the extension ID in its existing
  `extensions` configuration array.

Discovery does not mean enablement. Installing an extension wheel does not change
MDV behavior until its ID is configured.

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

### Startup behavior

| Installed | Configured | Result |
| --- | --- | --- |
| No | No | MDV starts normally. |
| Yes | No | MDV starts normally and does not import the extension. |
| Yes | Yes | MDV loads and activates the extension. |
| No | Yes | MDV stops with an `ExtensionError`. |

An enabled duplicate, malformed provider or startup failure also stops startup.
This fail-fast policy applies only when a deployment requested the extension;
normal deployments are unaffected by packages and IDs they did not enable.

## Application lifecycle

For each Flask application, MDV performs these steps:

1. Read and validate the configured extension IDs.
2. Discover installed `mdv.extensions` providers without loading their factories.
3. Resolve each configured ID from the built-in map or installed providers.
4. Construct each enabled extension once and validate optional capabilities.
5. Store the ordered mapping in `app.extensions`.
6. Register global routes and the generic navigation endpoint once.
7. Reuse the same objects when building project-specific server options.

This application ownership matters because Flask Blueprint registration is a
one-time side effect and an extension may own services or caches shared across
its routes.

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

The entry point should return a new object when its factory is called. MDV calls
the factory once for each Flask application; the extension should not implement
its own process-global singleton.

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

Hiding a catalog link is only a presentation decision. It does not protect the
Blueprint, its API routes or its static frontend. Route authorization remains
the responsibility of the host authentication contract and the extension's
operation-level checks.

## Backward compatibility

Existing built-in extensions continue to use the same configuration IDs and
three required methods. Installed discovery only handles configured IDs that are
not supplied by the built-in map. Navigation is optional, so extensions written
before navigation support remain valid.

This change does not alter MDV's database schema, project schema or extension
configuration shape. A deployment using only its previous built-in extension
array receives the same built-in extension objects in the same configured order,
with the existing always-on `ucsc_proxy` behavior retained.

## Responsibility boundaries

MDV owns:

- extension discovery, enablement, ordering and application-scoped reuse;
- validation of IDs, factories, the required interface and navigation metadata;
- global navigation rendering and MDV authentication identity;
- clear startup failure when a requested extension cannot be activated.

The external distribution owns:

- its Blueprint, frontend, business logic and third-party dependencies;
- authorization and validation for each operation;
- any extension-specific persistence and schema migration strategy;
- compatibility testing against the MDV revisions it supports.

The deployment owns:

- installing compatible external distributions;
- explicitly enabling the required IDs;
- pinning and testing the combined set of versions before delivery.
