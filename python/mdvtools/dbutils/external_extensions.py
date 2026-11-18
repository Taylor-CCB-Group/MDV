from __future__ import annotations

from collections.abc import Callable, Iterable
from importlib import metadata
from typing import Any

from mdvtools.server_extension import ExtensionError, MDVProjectServerExtension


EXTENSION_ENTRY_POINT_GROUP = "mdv.extensions"


ExternalExtensionProviders = dict[str, tuple[metadata.EntryPoint, ...]]


def discover_external_extensions() -> ExternalExtensionProviders:
    """Return installed MDV extension entry points grouped by extension ID."""

    discovered: dict[str, list[metadata.EntryPoint]] = {}
    for entry_point in _entry_points_for_group():
        discovered.setdefault(entry_point.name, []).append(entry_point)
    return {name: tuple(entry_points) for name, entry_points in discovered.items()}


def reject_external_provider(
    extension_id: str,
    providers: ExternalExtensionProviders,
) -> None:
    """Prevent an installed extension from shadowing an enabled built-in."""

    installed = providers.get(extension_id, ())
    if installed:
        provider_names = ", ".join(
            sorted(entry_point.value for entry_point in installed)
        )
        raise ExtensionError(
            f"Installed extension '{extension_id}' conflicts with an enabled built-in "
            f"extension: {provider_names}. Extension IDs must be unique."
        )


def load_external_extension(
    extension_id: str,
    providers: ExternalExtensionProviders,
) -> MDVProjectServerExtension:
    """Load and validate one explicitly enabled external extension."""

    installed = providers.get(extension_id, ())
    if not installed:
        raise ExtensionError(
            f"Extension '{extension_id}' is enabled but is neither built in "
            f"nor installed in entry-point group '{EXTENSION_ENTRY_POINT_GROUP}'."
        )
    if len(installed) > 1:
        provider_names = ", ".join(
            sorted(entry_point.value for entry_point in installed)
        )
        raise ExtensionError(
            f"Multiple installed distributions provide enabled extension "
            f"'{extension_id}': {provider_names}. Extension IDs must be unique."
        )

    entry_point = installed[0]
    try:
        factory = entry_point.load()
    except Exception as exc:
        raise ExtensionError(
            f"Extension '{extension_id}' failed while loading entry point "
            f"'{entry_point.value}'."
        ) from exc

    return _create_and_validate_extension(
        extension_id,
        factory,
        f"entry point '{entry_point.value}'",
    )


def _create_and_validate_extension(
    extension_id: str,
    factory: Callable[[], Any],
    source: str,
) -> MDVProjectServerExtension:
    if not callable(factory):
        raise ExtensionError(
            f"Extension '{extension_id}' from {source} must expose a zero-argument callable."
        )

    try:
        instance = factory()
    except Exception as exc:
        raise ExtensionError(
            f"Extension '{extension_id}' failed while being created from {source}."
        ) from exc

    missing_methods = [
        method_name
        for method_name in (
            "register_global_routes",
            "register_routes",
            "mutate_state_json",
        )
        if not callable(getattr(instance, method_name, None))
    ]
    if missing_methods:
        raise ExtensionError(
            f"Extension '{extension_id}' from {source} does not conform to "
            "MDVProjectServerExtension; missing callable methods: "
            + ", ".join(missing_methods)
        )

    return instance


def _entry_points_for_group() -> Iterable[metadata.EntryPoint]:
    return metadata.entry_points().select(group=EXTENSION_ENTRY_POINT_GROUP)
