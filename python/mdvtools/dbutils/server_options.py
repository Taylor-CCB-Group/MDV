from __future__ import annotations

from collections.abc import Iterable

from flask import Flask

from mdvtools.dbutils.extension_navigation import (
    register_extension_navigation_route,
    validate_extension_navigation,
)
from mdvtools.dbutils.external_extensions import (
    discover_external_extensions,
    load_external_extension,
    reject_external_provider,
)
from mdvtools.dbutils.project_manager_extension import ProjectManagerExtension
from mdvtools.llm.chat_server_extension import MDVProjectChatServerExtension
from mdvtools.server_extension import (
    ExtensionError,
    MDVProjectServerExtension,
    MDVServerOptions,
)
from mdvtools.ucsc_proxy_extension import UcscProxyServerExtension


ACTIVE_EXTENSIONS_KEY = "mdv.active_extensions"
GLOBAL_ROUTES_REGISTERED_KEY = "mdv.extension_global_routes_registered"


extension_classes: dict[str, type[MDVProjectServerExtension]] = {
    "chat": MDVProjectChatServerExtension,
    "project_manager": ProjectManagerExtension,
    # app-wide integrations that are safe to enable by default
    "ucsc_proxy": UcscProxyServerExtension,
}


def get_active_extensions(
    app: Flask,
) -> dict[str, MDVProjectServerExtension]:
    """Create configured extensions once and reuse them for this Flask app."""

    cached = app.extensions.get(ACTIVE_EXTENSIONS_KEY)
    if cached is not None:
        return cached

    configured_ids = _validate_configured_ids(app.config.get("extensions", []))
    external_providers = discover_external_extensions()

    # Preserve the existing always-on UCSC proxy behaviour.
    reject_external_provider("ucsc_proxy", external_providers)
    active: dict[str, MDVProjectServerExtension] = {
        "ucsc_proxy": UcscProxyServerExtension()
    }

    for extension_id in configured_ids:
        if extension_id == "ucsc_proxy":
            continue

        builtin_factory = extension_classes.get(extension_id)
        if builtin_factory is not None:
            reject_external_provider(extension_id, external_providers)
            instance = builtin_factory()
        else:
            instance = load_external_extension(extension_id, external_providers)

        validate_extension_navigation(extension_id, instance)
        active[extension_id] = instance

    app.extensions[ACTIVE_EXTENSIONS_KEY] = active
    return active


def get_active_extension(
    app: Flask, extension_id: str
) -> MDVProjectServerExtension | None:
    """Return one active app-local extension instance by configured ID."""

    return get_active_extensions(app).get(extension_id)


def register_global_routes_for_extensions(
    app: Flask,
) -> dict[str, MDVProjectServerExtension]:
    """Register host navigation and extension global routes exactly once per app."""

    active = get_active_extensions(app)
    if app.extensions.get(GLOBAL_ROUTES_REGISTERED_KEY):
        return active

    register_extension_navigation_route(app, active)
    for extension_id, extension in active.items():
        try:
            extension.register_global_routes(app, app.config)
        except Exception as exc:
            raise ExtensionError(
                f"Extension '{extension_id}' failed while registering global routes."
            ) from exc

    app.extensions[GLOBAL_ROUTES_REGISTERED_KEY] = True
    return active


def get_server_options_for_db_projects(app: Flask) -> MDVServerOptions:
    """Return database-project options using the app's active extensions."""

    return MDVServerOptions(
        open_browser=False,
        backend_db=True,
        app=app,
        extensions=list(get_active_extensions(app).values()),
        websocket=True,
    )


def _validate_configured_ids(configured_ids: object) -> tuple[str, ...]:
    if not isinstance(configured_ids, list):
        raise ExtensionError("MDV config 'extensions' must be a list of extension IDs.")
    if any(
        not isinstance(extension_id, str) or not extension_id
        for extension_id in configured_ids
    ):
        raise ExtensionError(
            "MDV config 'extensions' must contain non-empty string IDs."
        )

    duplicates = _duplicates(configured_ids)
    if duplicates:
        raise ExtensionError(
            "MDV config contains duplicate extension IDs: " + ", ".join(duplicates)
        )
    return tuple(configured_ids)


def _duplicates(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return sorted(duplicates)
