import json

from mdvtools.llm import chat_cli
from mdvtools.mdvproject import MDVProject


def _make_project(tmp_path):
    project_dir = tmp_path / "proj"
    project = MDVProject(str(project_dir))
    project.datasources = [{"name": "cells", "columns": [{"field": "cell_id"}]}]
    return project


def test_run_chat_once_success_with_output_dir(tmp_path, monkeypatch):
    project = _make_project(tmp_path)

    class FakeProjectChat:
        def __init__(self, incoming_project):
            self.project = incoming_project
            self.init_error = False

        def ask_question(self, _chat_request):
            self.project.set_view("v_test", {"initialCharts": {"cells": []}})
            return {
                "code": "I ran code\n```python\nprint('ok')\n```\n",
                "view_name": "v_test",
                "error": False,
                "message": "Success",
                "verification": None,
                "data_preview": None,
                "needs_refresh": False,
            }

    monkeypatch.setattr(chat_cli, "ProjectChat", FakeProjectChat)
    output_dir = tmp_path / "debug"
    result, exit_code = chat_cli.run_chat_once(
        project_path=str(project.dir),
        prompt="build a view",
        output_dir=str(output_dir),
    )

    assert exit_code == 0
    assert result["success"] is True
    assert result["view_name"] == "v_test"
    assert result["debug_output_dir"] == str(output_dir.resolve())

    generated_code = (output_dir / "generated_code.py").read_text(encoding="utf-8")
    assert "print('ok')" in generated_code
    payload = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    assert payload["result"]["success"] is True
    view_payload = json.loads((output_dir / "view_snapshot.json").read_text(encoding="utf-8"))
    assert view_payload["view_name"] == "v_test"
    assert view_payload["view"]["initialCharts"] == {"cells": []}


def test_run_chat_once_failure_with_debug_output(tmp_path):
    missing_project = tmp_path / "does-not-exist"
    output_dir = tmp_path / "debug-failure"

    result, exit_code = chat_cli.run_chat_once(
        project_path=str(missing_project),
        prompt="build a view",
        output_dir=str(output_dir),
    )

    assert exit_code == 2
    assert result["success"] is False
    assert isinstance(result["duration_seconds"], float)
    assert result["duration_seconds"] >= 0
    assert result["view_name"] is None
    assert result["debug_output_dir"] == str(output_dir.resolve())
    assert (output_dir / "result.json").exists()


def test_main_json_output_schema(tmp_path, monkeypatch, capsys):
    project = _make_project(tmp_path)

    class FakeProjectChat:
        def __init__(self, incoming_project):
            self.project = incoming_project
            self.init_error = False

        def ask_question(self, _chat_request):
            return {
                "code": "```python\nprint('ok')\n```",
                "view_name": "v_json",
                "error": False,
                "message": "Success",
                "verification": None,
                "data_preview": None,
                "needs_refresh": False,
            }

    monkeypatch.setattr(chat_cli, "ProjectChat", FakeProjectChat)

    exit_code = chat_cli.main(
        ["--project", str(project.dir), "--prompt", "run", "--json"]
    )
    stdout = capsys.readouterr().out
    payload = json.loads(stdout)

    assert exit_code == 0
    assert payload["success"] is True
    assert payload["project_path"] == str(project.dir)
    assert payload["views_file"].endswith("views.json")
    assert payload["view_name"] == "v_json"
    assert "debug_output_dir" in payload


def test_run_chat_once_installs_socketio_shim_when_missing(tmp_path, monkeypatch):
    project = _make_project(tmp_path)
    from mdvtools import websocket as mdv_websocket
    monkeypatch.setattr(mdv_websocket, "socketio", None)

    class FakeProjectChat:
        def __init__(self, incoming_project):
            self.project = incoming_project
            self.init_error = False

        def ask_question(self, _chat_request):
            return {
                "code": "```python\nprint('ok')\n```",
                "view_name": "shim_view",
                "error": False,
                "message": "Success",
                "verification": None,
                "data_preview": None,
                "needs_refresh": False,
            }

    monkeypatch.setattr(chat_cli, "ProjectChat", FakeProjectChat)
    result, exit_code = chat_cli.run_chat_once(
        project_path=str(project.dir),
        prompt="run",
    )

    assert exit_code == 0
    assert result["success"] is True
    assert mdv_websocket.socketio is not None


def test_run_chat_once_passes_model_id_to_ask_question(tmp_path, monkeypatch):
    project = _make_project(tmp_path)
    captured: dict[str, object] = {}

    class FakeProjectChat:
        def __init__(self, incoming_project):
            self.project = incoming_project
            self.init_error = False

        def ask_question(self, chat_request):
            captured["chat_request"] = chat_request
            return {
                "code": "```python\nprint('ok')\n```",
                "view_name": "v_model",
                "error": False,
                "message": "Success",
                "verification": None,
                "data_preview": None,
                "needs_refresh": False,
            }

    monkeypatch.setattr(chat_cli, "ProjectChat", FakeProjectChat)
    monkeypatch.setattr(
        chat_cli,
        "resolve_chat_model_id",
        lambda model: "openai:chat:gpt-4o-mini" if model == "gpt-4o-mini" else None,
    )

    result, exit_code = chat_cli.run_chat_once(
        project_path=str(project.dir),
        prompt="run",
        model_id="gpt-4o-mini",
    )

    assert exit_code == 0
    assert captured["chat_request"]["model_id"] == "openai:chat:gpt-4o-mini"
    assert result["model_id"] == "openai:chat:gpt-4o-mini"


def test_main_rejects_unknown_model(tmp_path, monkeypatch, capsys):
    project = _make_project(tmp_path)
    monkeypatch.setattr(
        chat_cli,
        "resolve_chat_model_id",
        lambda _model: (_ for _ in ()).throw(ValueError("Unknown chat model 'bad'")),
    )

    exit_code = chat_cli.main(
        [
            "--project",
            str(project.dir),
            "--prompt",
            "run",
            "--model",
            "bad",
        ]
    )
    assert exit_code == 2
    assert "Unknown chat model" in capsys.readouterr().err


def test_run_chat_once_renames_view_with_view_name(tmp_path, monkeypatch):
    project = _make_project(tmp_path)

    class FakeProjectChat:
        def __init__(self, incoming_project):
            self.project = incoming_project
            self.init_error = False

        def ask_question(self, _chat_request):
            self.project.set_view("auto_view", {"initialCharts": {"cells": []}})
            return {
                "code": "```python\nprint('ok')\n```",
                "view_name": "auto_view",
                "error": False,
                "message": "Success",
                "verification": None,
                "data_preview": None,
                "needs_refresh": False,
            }

    monkeypatch.setattr(chat_cli, "ProjectChat", FakeProjectChat)
    result, exit_code = chat_cli.run_chat_once(
        project_path=str(project.dir),
        prompt="run",
        view_name="custom_view",
    )

    assert exit_code == 0
    assert result["success"] is True
    assert result["view_name"] == "custom_view"
    assert project.get_view("custom_view") is not None
    assert project.get_view("auto_view") is None
