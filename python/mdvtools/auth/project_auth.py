from flask import jsonify, session, current_app
from datetime import datetime, timedelta
from mdvtools.auth.authutils import user_project_cache
from mdvtools.dbutils.dbservice import ProjectService
from mdvtools.logging_config import get_logger

logger = get_logger(__name__)

TIMESTAMP_UPDATE_INTERVAL = timedelta(hours=1)

def project_pre_dispatch_checks(project_id: str, options: dict):
    """
    Performs checks before dispatching a request to a project route.
    This function is intended to be used as a 'before_request' hook.

    - Checks if the project exists.
    - Updates the project's 'accessed_timestamp'.
    - Performs authentication and authorization checks based on the route's options.

    :param project_id: The ID of the project being accessed.
    :param options: The options dictionary for the matched route (e.g., {'access_level': 'editable'}).
    :return: A Flask Response object to short-circuit the request in case of an error, or None to continue.
    """
    project = ProjectService.get_project_by_id(project_id)
    if not project:
        logger.error(f"In pre-dispatch check: Project with ID {project_id} not found.")
        return jsonify({"error": f"Project with ID {project_id} not found"}), 404

    # Update the accessed timestamp only if the last update was more than TIMESTAMP_UPDATE_INTERVAL ago
    if not project.accessed_timestamp or (datetime.now() - project.accessed_timestamp > TIMESTAMP_UPDATE_INTERVAL):
        try:
            ProjectService.set_project_accessed_timestamp(project_id)
        except Exception as e:
            logger.exception(f"Error updating accessed timestamp for project {project_id}: {e}")
            return jsonify({"error": "Failed to update project accessed timestamp"}), 500

    # If auth is not enabled for the app, we can skip the user-specific checks
    if not current_app.config.get("ENABLE_AUTH", False):
        if options.get('access_level') == 'editable' and (
            not project.access_level or project.access_level != 'editable'
        ):
            return jsonify({"error": "This project is read-only and cannot be modified."}), 403
        return None # No more checks needed

    # --- Auth-enabled checks below ---
    user = session.get('user')
    user_id = user.get("id") if user else None
    if not user_id:
        return jsonify({"error": "Authentication required."}), 401

    try:
        project_permissions = user_project_cache.get(user_id, {}).get(int(project_id))
    except ValueError as e:
        # project_id is passed as a string, conversion to int may fail
        logger.exception(f"Error getting project permissions for user {user_id} and project {project_id}: {e}")
        return jsonify({"error": "Failed to get project permissions"}), 500

    # Check for access level only if specified in options
    if options.get('access_level') == 'editable':
        # Check permission from the database model
        if project.access_level != 'editable':
            return jsonify({"error": "This project is read-only and cannot be modified."}), 403
        
        # Check user-specific write permissions
        if not project_permissions or not (project_permissions.get("is_owner") or project_permissions.get("can_write")):
            return jsonify({"error": "User does not have the required write permissions for this project."}), 403

        # On successful write access, update the project's update timestamp
        ProjectService.set_project_update_timestamp(project_id)

    return None # All checks passed 