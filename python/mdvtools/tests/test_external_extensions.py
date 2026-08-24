from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pytest
from flask import Flask

from mdvtools.dbutils import external_extensions
from mdvtools.dbutils.project_manager_extension import ProjectManagerExtension
from mdvtools.dbutils.server_options import (
    ACTIVE_EXTENSIONS_KEY,
    get_active_extensions,
    get_server_options_for_db_projects,
    register_global_routes_for_extensions,
)
from mdvtools.llm.chat_server_extension import MDVProjectChatServerExtension
from mdvtools.server_extension import ExtensionError, ExtensionNavigation
from mdvtools.ucsc_proxy_extension import UcscProxyServerExtension


@dataclass
class FakeEntryPoint:
    name: str
    value: str
    factory: Any
    load_error: Exception | None = None
    load_count: int = 0

    def load(self):
        self.load_count += 1
        if self.load_error is not None:
            raise self.load_error
        return self.factory


class FakeExtension:
    def __init__(
        self,
        extension_id: str = "example_extension",
        navigation: ExtensionNavigation | None = None,
        fail_registration: bool = False,
    ) -> None:
        self.extension_id = extension_id
        self.navigation = navigation
        self.fail_registration = fail_registration
        self.global_registration_count = 0

    def register_global_routes(self, app: Flask, config: dict) -> None:
        self.global_registration_count += 1
        if self.fail_registration:
            raise ValueError("test registration failed")
        app.add_url_rule(
            f"/{self.extension_id}/status",
            endpoint=f"{self.extension_id}_status",
            view_func=lambda: {"extension": self.extension_id, "status": "ok"},
        )

    def register_routes(self, project, project_bp) -> None:
        return None

    def mutate_state_json(self, state_json: dict, project, app: Flask) -> None:
        return None


def use_entry_points(monkeypatch, *entry_points: FakeEntryPoint) -> None:
    monkeypatch.setattr(
        external_extensions,
        "_entry_points_for_group",
        lambda: entry_points,
    )


def extension_factory(
    extension_id: str = "example_extension",
    navigation: Any = None,
    fail_registration: bool = False,
):
    return lambda: FakeExtension(extension_id, navigation, fail_registration)


def test_only_enabled_external_extensions_are_loaded_and_ordered(monkeypatch) -> None:
    enabled = FakeEntryPoint(
        "example_extension",
        "example:create_extension",
        extension_factory(
            navigation=ExtensionNavigation(
                label="Example Extension",
                url="/example-extension/",
            )
        ),
    )
    disabled = FakeEntryPoint(
        "installed_only",
        "other:create_extension",
        lambda: (_ for _ in ()).throw(ValueError("must stay unloaded")),
    )
    use_entry_points(monkeypatch, disabled, enabled)
    app = Flask(__name__)
    app.secret_key = "test"
    app.config["extensions"] = ["project_manager", "example_extension"]

    first = get_active_extensions(app)
    second = get_active_extensions(app)
    first_options = get_server_options_for_db_projects(app)
    second_options = get_server_options_for_db_projects(app)
    register_global_routes_for_extensions(app)
    register_global_routes_for_extensions(app)
    example_instance = first["example_extension"]

    assert list(first) == [
        "ucsc_proxy",
        "project_manager",
        "example_extension",
    ]
    assert first is second
    assert first_options.extensions == second_options.extensions
    assert first_options.extensions[-1] is example_instance
    assert app.extensions[ACTIVE_EXTENSIONS_KEY] is first
    assert enabled.load_count == 1
    assert disabled.load_count == 0
    assert isinstance(example_instance, FakeExtension)
    assert example_instance.global_registration_count == 1
    assert app.test_client().get("/example_extension/status").get_json() == {
        "extension": "example_extension",
        "status": "ok",
    }
    assert app.test_client().get("/extension_navigation").get_json() == {
        "extensions": [
            {
                "id": "example_extension",
                "label": "Example Extension",
                "url": "/example-extension/",
            }
        ]
    }


def test_duplicate_disabled_entry_points_remain_unloaded(monkeypatch) -> None:
    first = FakeEntryPoint(
        "unused",
        "one:create_extension",
        extension_factory("unused"),
    )
    second = FakeEntryPoint(
        "unused",
        "two:create_extension",
        extension_factory("unused"),
    )
    use_entry_points(monkeypatch, first, second)
    app = Flask(__name__)
    app.config["extensions"] = []

    active = get_active_extensions(app)

    assert list(active) == ["ucsc_proxy"]
    assert first.load_count == 0
    assert second.load_count == 0


def test_existing_builtin_extensions_remain_ordered_and_reused(monkeypatch) -> None:
    use_entry_points(monkeypatch)
    app = Flask(__name__)
    app.config["extensions"] = ["chat", "project_manager", "ucsc_proxy"]

    active = get_active_extensions(app)
    options = get_server_options_for_db_projects(app)

    assert list(active) == ["ucsc_proxy", "chat", "project_manager"]
    assert isinstance(active["ucsc_proxy"], UcscProxyServerExtension)
    assert isinstance(active["chat"], MDVProjectChatServerExtension)
    assert isinstance(active["project_manager"], ProjectManagerExtension)
    assert options.extensions == list(active.values())
    assert get_active_extensions(app) is active


def test_missing_enabled_extension_is_rejected(monkeypatch) -> None:
    use_entry_points(monkeypatch)
    app = Flask(__name__)
    app.config["extensions"] = ["missing_extension"]

    with pytest.raises(
        ExtensionError, match="enabled but is neither built in nor installed"
    ):
        get_active_extensions(app)


def test_duplicate_installed_extension_ids_are_rejected_when_enabled(
    monkeypatch,
) -> None:
    first = FakeEntryPoint(
        "example_extension", "one:create_extension", extension_factory()
    )
    second = FakeEntryPoint(
        "example_extension", "two:create_extension", extension_factory()
    )
    use_entry_points(monkeypatch, first, second)
    app = Flask(__name__)
    app.config["extensions"] = ["example_extension"]

    with pytest.raises(ExtensionError, match="Multiple installed distributions"):
        get_active_extensions(app)

    assert first.load_count == 0
    assert second.load_count == 0


def test_duplicate_configured_extension_ids_are_rejected(monkeypatch) -> None:
    use_entry_points(
        monkeypatch,
        FakeEntryPoint(
            "example_extension", "example:create_extension", extension_factory()
        ),
    )
    app = Flask(__name__)
    app.config["extensions"] = ["example_extension", "example_extension"]

    with pytest.raises(ExtensionError, match="duplicate extension IDs"):
        get_active_extensions(app)


def test_external_extension_cannot_shadow_enabled_builtin(monkeypatch) -> None:
    entry_point = FakeEntryPoint(
        "project_manager",
        "private:create_extension",
        extension_factory("project_manager"),
    )
    use_entry_points(monkeypatch, entry_point)
    app = Flask(__name__)
    app.config["extensions"] = ["project_manager"]

    with pytest.raises(ExtensionError, match="conflicts with an enabled built-in"):
        get_active_extensions(app)

    assert entry_point.load_count == 0


@pytest.mark.parametrize(
    "entry_point, message",
    [
        (
            FakeEntryPoint(
                "example_extension",
                "example:create_extension",
                extension_factory(),
                load_error=ValueError("load failed"),
            ),
            "failed while loading entry point",
        ),
        (
            FakeEntryPoint(
                "example_extension",
                "example:not_callable",
                object(),
            ),
            "zero-argument callable",
        ),
        (
            FakeEntryPoint(
                "example_extension",
                "example:create_extension",
                lambda: object(),
            ),
            "does not conform",
        ),
    ],
)
def test_malformed_enabled_extensions_are_rejected(
    monkeypatch, entry_point, message
) -> None:
    use_entry_points(monkeypatch, entry_point)
    app = Flask(__name__)
    app.config["extensions"] = ["example_extension"]

    with pytest.raises(ExtensionError, match=message):
        get_active_extensions(app)


def test_broken_global_registration_is_rejected(monkeypatch) -> None:
    use_entry_points(
        monkeypatch,
        FakeEntryPoint(
            "example_extension",
            "example:create_extension",
            extension_factory(fail_registration=True),
        ),
    )
    app = Flask(__name__)
    app.config["extensions"] = ["example_extension"]

    with pytest.raises(ExtensionError, match="failed while registering global routes"):
        register_global_routes_for_extensions(app)


@pytest.mark.parametrize(
    "navigation, message",
    [
        (ExtensionNavigation(label="", url="/example/"), "label must not be empty"),
        (
            ExtensionNavigation(label="Example", url="https://example.com"),
            "app-relative path",
        ),
        (
            ExtensionNavigation(label="Example", url="//example.com/path"),
            "app-relative path",
        ),
        ({"label": "Example", "url": "/example/"}, "must be ExtensionNavigation"),
    ],
)
def test_invalid_navigation_is_rejected(monkeypatch, navigation, message) -> None:
    use_entry_points(
        monkeypatch,
        FakeEntryPoint(
            "example_extension",
            "example:create_extension",
            extension_factory(navigation=navigation),
        ),
    )
    app = Flask(__name__)
    app.config["extensions"] = ["example_extension"]

    with pytest.raises(ExtensionError, match=message):
        get_active_extensions(app)


def test_admin_navigation_is_filtered_by_mdv_identity(monkeypatch) -> None:
    use_entry_points(
        monkeypatch,
        FakeEntryPoint(
            "admin_extension",
            "admin:create_extension",
            extension_factory(
                "admin_extension",
                ExtensionNavigation(
                    label="Admin",
                    url="/admin/",
                    requires_admin=True,
                ),
            ),
        ),
    )
    app = Flask(__name__)
    app.secret_key = "test"
    app.config.update(ENABLE_AUTH=True, extensions=["admin_extension"])
    register_global_routes_for_extensions(app)
    client = app.test_client()

    assert client.get("/extension_navigation").get_json() == {"extensions": []}

    with client.session_transaction() as flask_session:
        flask_session["user"] = {"id": 1, "is_admin": False}
    assert client.get("/extension_navigation").get_json() == {"extensions": []}

    with client.session_transaction() as flask_session:
        flask_session["user"] = {"id": 2, "is_admin": True}
    assert client.get("/extension_navigation").get_json() == {
        "extensions": [{"id": "admin_extension", "label": "Admin", "url": "/admin/"}]
    }
