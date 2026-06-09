import os

from flask import session

from mdvtools.auth import authutils
from mdvtools.mdvproject import MDVProject
from mdvtools.dbutils.dbservice import ProjectService
from mdvtools.project_router import ProjectBlueprint


def refresh_auth_cache(enable_auth: bool):
    if enable_auth:
        authutils.cache_user_projects()


def get_current_user_id(enable_auth: bool):
    if not enable_auth:
        return None
    user = session.get("user")
    if not user:
        raise PermissionError("User not found in session.")
    return user["id"]


def validate_owned_projects(project_ids: list[int], enable_auth: bool):
    if not enable_auth:
        return
    user_id = get_current_user_id(enable_auth)
    user_projects = authutils.user_project_cache.get(user_id, {})
    unauthorized_ids = [
        project_id
        for project_id in project_ids
        if not user_projects.get(project_id, {}).get("is_owner", False)
    ]
    if unauthorized_ids:
        raise PermissionError(f"Only the project owner can delete projects: {unauthorized_ids}")


def soft_delete_projects(project_ids: list[int], enable_auth: bool):
    validate_owned_projects(project_ids, enable_auth)
    deleted_projects = ProjectService.soft_delete_projects(project_ids)
    for project in deleted_projects:
        ProjectBlueprint.blueprints.pop(str(project.id), None)
    refresh_auth_cache(enable_auth)
    return deleted_projects


def restore_deleted_project(project_id: int, enable_auth: bool, app):
    validate_owned_projects([project_id], enable_auth)
    project = ProjectService.get_project_by_id(project_id)
    if not project:
        raise ValueError(f"Project with ID {project_id} not found.")
    if not project.is_deleted:
        raise ValueError(f"Project with ID {project_id} is not in the recycle bin.")
    if not project.path or not os.path.exists(project.path):
        raise ValueError(f"Project with ID {project_id} is missing its filesystem path.")

    try:
        mdv_project = MDVProject(dir=project.path, id=str(project.id), backend_db=True)
        mdv_project.serve(app=app, open_browser=False, backend_db=True)
        restored_project = ProjectService.restore_deleted_project(project_id)
        refresh_auth_cache(enable_auth)
        return restored_project
    except Exception:
        ProjectBlueprint.blueprints.pop(str(project_id), None)
        raise
