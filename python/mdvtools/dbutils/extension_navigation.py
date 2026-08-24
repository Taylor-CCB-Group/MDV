from __future__ import annotations

from collections.abc import Mapping
from typing import Any
from urllib.parse import urlsplit

from flask import Flask, session

from mdvtools.server_extension import (
    ExtensionError,
    ExtensionNavigation,
    MDVProjectServerExtension,
)


def validate_extension_navigation(extension_id: str, instance: Any) -> None:
    """Validate optional catalog navigation supplied by an extension."""

    navigation = getattr(instance, "navigation", None)
    if navigation is None:
        return
    if not isinstance(navigation, ExtensionNavigation):
        raise ExtensionError(
            f"Extension '{extension_id}' navigation must be ExtensionNavigation."
        )
    if not navigation.label.strip():
        raise ExtensionError(
            f"Extension '{extension_id}' navigation label must not be empty."
        )
    parsed_url = urlsplit(navigation.url)
    if not navigation.url.startswith("/") or parsed_url.scheme or parsed_url.netloc:
        raise ExtensionError(
            f"Extension '{extension_id}' navigation URL must be an app-relative path."
        )
    if not isinstance(navigation.requires_admin, bool):
        raise ExtensionError(
            f"Extension '{extension_id}' navigation requires_admin must be a boolean."
        )


def register_extension_navigation_route(
    app: Flask,
    active_extensions: Mapping[str, MDVProjectServerExtension],
) -> None:
    """Expose navigation metadata from the app's active extensions."""

    def extension_navigation():
        auth_enabled = app.config.get("ENABLE_AUTH", False)
        user = session.get("user") if auth_enabled else None
        is_admin = not auth_enabled or bool(user and user.get("is_admin"))

        entries = []
        for extension_id, extension in active_extensions.items():
            navigation = getattr(extension, "navigation", None)
            if navigation is None:
                continue
            if navigation.requires_admin and not is_admin:
                continue
            entries.append(
                {
                    "id": extension_id,
                    "label": navigation.label,
                    "url": navigation.url,
                }
            )
        return {"extensions": entries}

    app.add_url_rule(
        "/extension_navigation",
        endpoint="mdv_extension_navigation",
        view_func=extension_navigation,
    )
