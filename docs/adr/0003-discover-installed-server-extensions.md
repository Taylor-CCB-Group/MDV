# Discover installed packages through the existing server extension interface

**Status:** proposed

## Context

MDV already has a server extension interface used by built-in features. An
extension can register routes once for a Flask application, add project routes,
and mutate project state returned to the browser. Built-in extension classes were
selected from the deployment's existing `extensions` array, but the available
classes were hardcoded inside MDV.

That hardcoding prevents a separately distributed package from participating in
the lifecycle unless public MDV imports the package directly. This is unsuitable
for optional or private products: public MDV must remain installable without
them, and their source and dependencies must remain in their own distributions.

## Decision

MDV will evolve the existing server extension system instead of introducing a
parallel plugin system.

- The existing `extensions` configuration array remains the explicit list of
  enabled extension IDs.
- Built-in extensions remain in MDV's built-in map.
- Separately installed Python distributions advertise zero-argument factories
  through the `mdv.extensions` entry-point group.
- MDV combines built-in and installed providers while resolving the configured
  IDs. It does not import an installed provider unless that ID is enabled.
- Each enabled extension is instantiated once per Flask application. MDV stores
  the resolved mapping in `app.extensions` and reuses the instances for global
  and project behavior.
- Catalog navigation is an optional capability on the existing interface. It is
  not a new required method.
- Configuration and provider errors are detected during application startup.

## Why the existing interface is the right boundary

The required lifecycle already exists. A global-only product can register a
normal Flask Blueprint in `register_global_routes()` and implement the two
project methods as no-ops. Adding another plugin lifecycle would create two ways
to configure, initialise and reason about the same kind of server integration.

Reusing the interface also preserves built-in behavior. Existing extensions keep
their classes and methods, while installed packages gain a packaging-level
discovery mechanism. No database schema or project-state schema is introduced by
this decision.

## Why Python entry points are used

Entry points are metadata in an installed Python distribution. They allow a
package to declare a provider without MDV knowing the package's import path or
source location.

This is preferable to:

- hardcoding optional packages in public MDV, which creates a public dependency
  on packages that may be absent or private;
- putting Python import strings in deployment JSON, which exposes implementation
  paths and postpones validation until import time;
- scanning directories or source repositories, which couples runtime discovery
  to a particular checkout layout;
- automatically enabling every installed provider, which makes installing a
  wheel change application behavior without an explicit configuration decision.

The entry-point name is the extension ID. IDs must be unique so configuration,
errors, navigation and logs all refer to one unambiguous provider.

## Installation and enablement are separate

An installed distribution makes an extension *available*. The `extensions`
array makes it *enabled* for one deployment.

This distinction gives normal public deployments and premium deployments the
same MDV code with different installed packages and configuration:

| Situation | Startup behavior | Reason |
| --- | --- | --- |
| Package is not installed and ID is not configured | Start normally | The deployment did not request the extension. |
| Package is installed and ID is not configured | Start normally without importing it | Installation alone must not activate behavior. |
| Package is installed and ID is configured | Load and activate it | Both availability and operator intent are present. |
| ID is configured but no provider exists | Stop with a clear error | A requested capability is missing or misspelled. |
| Enabled provider is duplicated, malformed or fails startup | Stop with a clear error | Continuing would leave a partially configured application. |

Silently ignoring a configured but missing extension was rejected because a
misspelled or incomplete premium deployment could appear healthy while omitting
required routes and controls.

## Why instances are application-scoped

Global Flask routes can only be registered safely once, and extensions may own
service objects or caches that must be shared by their global and project
behavior. Constructing fresh instances at different call sites would make their
state inconsistent and could repeat side effects.

The Flask application is therefore the ownership boundary. A second Flask
application in the same process receives its own instances, while every use
inside one application receives the same objects.

## Why navigation is optional

Existing server extensions predate catalog navigation. Making navigation a new
required method would break otherwise valid extensions. Instead, an extension
may expose validated `ExtensionNavigation` metadata. The host renders that
metadata generically and can omit administrator links for users who are not MDV
administrators.

Navigation visibility is not authorization. Extension routes must still use
MDV's authenticated identity and authorization helpers.

## Considered and rejected

- **A separate plugin configuration and lifecycle.** Duplicates the existing
  extension seam and leaves built-in and installed integrations behaving
  differently.
- **Hardcode the external package in MDV.** Makes an optional/private package a
  public MDV concern and prevents MDV from starting when it is absent.
- **Silently ignore configured missing extensions.** Hides deployment errors.
- **Automatically activate installed extensions.** Installation would have an
  unexpected runtime effect and make deployment intent unclear.
- **Create an extension instance at every existing call site.** Risks repeated
  route registration and splits extension-owned state.
- **Make navigation part of the required protocol.** Breaks legacy extensions
  that do not need catalog UI.

## Consequences

- Public MDV knows only the generic interface and entry-point group; it does not
  import private product code.
- External packages can package and import their own dependencies normally.
- Existing built-in IDs and configuration continue to work.
- Global-only extensions use two explicit no-op project methods.
- Invalid enabled extensions fail early rather than producing a partially
  working deployment.
- The deployment artifact is responsible for installing compatible external
  wheels and enabling their IDs. Installation/version pinning remains a delivery
  concern rather than an extension-lifecycle concern.
