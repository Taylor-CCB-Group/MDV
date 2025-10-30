import os
from typing import Any
from mdvtools.logging_config import get_logger
from mdvtools.dbutils.project_manager_extension import ProjectManagerExtension
logger = get_logger(__name__)

# Module-level constants and functions that can be imported
REQUIRED_FILES = {"views.json", "state.json", "datasources.json"}

def get_project_thumbnail(project_path):
    """Extract the first available viewImage from a project's views."""
    try:
        from mdvtools.mdvproject import MDVProject
        mdv_project = MDVProject(project_path)
        return next((v["viewImage"] for v in mdv_project.views.values() if "viewImage" in v), None)
    except Exception as e:
        logger.exception(f"Error extracting thumbnail for project at {project_path}: {e}")
        return None

def find_root_prefix(names):
    # Check if all required files are in the root of archive
    if REQUIRED_FILES.issubset(set(os.path.basename(n) for n in names if "/" not in n)):
        return ""
    # Check one level below if required files exist
    dirs = {n.split("/", 1)[0] for n in names if "/" in n}
    for d in dirs:
        files_in_d = {os.path.basename(n) for n in names if n.startswith(f"{d}/")}
        if REQUIRED_FILES.issubset(files_in_d):
            return d + "/"
    return None

def register_routes(app, ENABLE_AUTH):
    from flask import abort, request, jsonify, session, redirect, url_for, render_template
    from mdvtools.auth.authutils import active_projects_cache, user_project_cache, user_cache, all_users_cache, cache_user_projects
    from mdvtools.dbutils.mdv_server_app import serve_projects_from_filesystem
    from mdvtools.dbutils.dbservice import ProjectService, UserProjectService
    import os
    import shutil
    from mdvtools.mdvproject import MDVProject
    from mdvtools.project_router import ProjectBlueprint
    from mdvtools.dbutils.dbmodels import User
    import tempfile
    
    """Register routes with the Flask app."""
    logger.info("Registering routes...")
    # Note: Project management routes (create, import, export, delete, rename, access, share)
    # are now handled by the ProjectManagerExtension

    try:
        @app.route('/')
        def index():
            try:
                return render_template('index.html')
            except Exception as e:
                logger.exception(f"Error rendering index: {e}")
                return jsonify({"error": str(e)}), 500

        logger.info("Route registered: /")

        @app.route('/login_dev')
        def login_dev():
            return render_template('login.html')
        logger.info("Route registered: /login_dev")

        # This is a relative route, we need to change it at some point
        @app.route('/api_root')
        def get_mdv_api_root():
            root = os.environ.get("MDV_API_ROOT", "/") or "/"
            if not root.startswith("/"):
                root = "/" + root
            if not root.endswith("/"):
                root = root + "/"
            return jsonify({"mdv_api_root": root})
        logger.info("Route registered: /api_root")

        @app.route('/enable_auth')
        def get_enable_auth():
            enable_auth = os.environ.get("ENABLE_AUTH") or "0"
            is_auth_enabled = True if enable_auth == "1" else False
            return jsonify({"enable_auth": is_auth_enabled})
        logger.info("Route registered: /enable_auth")

        @app.route('/extension_config', methods=['GET'])
        def extension_config():
            """Return a static-like map of extension name to its API flags.
            """
            try:
                from mdvtools.dbutils.project_manager_extension import ProjectManagerExtension

                # Determine if project_manager is enabled in config.json
                enabled_extensions = app.config.get('extensions', []) or []
                pm_enabled = 'project_manager' in enabled_extensions

                # Get the canonical set of keys from the extension, then override values
                pm = ProjectManagerExtension()
                ext_any: Any = pm
                pm_config_true = ext_any.get_session_config()  # type: ignore[attr-defined]

                if pm_enabled:
                    pm_config = pm_config_true
                else:
                    # Same keys, all set to False
                    pm_config = {k: False for k in pm_config_true.keys()}

                return jsonify({
                    "project_manager": pm_config
                })
            except Exception as e:
                logger.exception(f"Error in /extension_config: {e}")
                return jsonify({}), 500

        @app.route('/rescan_projects', methods=['GET'])
        def rescan_projects():
            if ENABLE_AUTH:
                user = session.get('user')
                if not user:
                    # If user is not found in session, raise an error or redirect.
                    raise ValueError("User not found in session.")
                
                is_admin = user.get("is_admin", False)
                if not is_admin:
                    # If user isn't admin, abort with 403 (Forbidden).
                    abort(403)  # Forbidden

            try:
                # Serve the projects after checking authentication and admin privileges
                created_ids = serve_projects_from_filesystem(app, app.config["projects_base_dir"])

                # If auth is enabled and there are new projects, grant owner to current admin and refresh cache
                if ENABLE_AUTH and created_ids:
                    try:
                        user_id = user.get("id") if user else None
                        if user_id is not None:
                            for pid in created_ids:
                                UserProjectService.add_or_update_user_project(user_id=user_id, project_id=pid, is_owner=True)
                        # Refresh caches so /projects reflects new permissions
                        cache_user_projects()
                    except Exception as perm_e:
                        logger.exception(f"Error assigning permissions for new projects {created_ids}: {perm_e}")
            except Exception as e:
                # Handle potential errors while serving the projects
                logger.exception(f"Error while serving the projects: {e}")
                abort(500, description="Error while serving the projects.")

            # If everything goes well, redirect to the 'index' page
            return redirect(url_for('index'))
        
        def get_project_owners(project_id):
            owners = []
            for user_id, projects in user_project_cache.items():
                perms = projects.get(project_id)
                if perms and perms.get("is_owner"):
                    user_data = next((u for u in all_users_cache if u["id"] == user_id), None)
                    if user_data:
                        owners.append(user_data["email"])
            return sorted(owners)

        @app.route('/projects')
        def get_projects():
            """
            Fetches the list of active projects from cache or database depending on the authentication.
            """
            logger.info(" /projects queried...")

            try:
                # Step 1: Get projects directly from cache
                if ENABLE_AUTH:
                    active_projects = active_projects_cache
                else:
                    # If authentication is disabled, get all projects from the database
                    active_projects = ProjectService.get_active_projects()

                if ENABLE_AUTH:
                    user = session.get('user')
                    if not user:
                        raise ValueError("User not found in session.")
                    # Filter by user permissions if authentication is enabled
                    user_id = user["id"]
                    user_projects = user_project_cache.get(user_id)
                    
                    # If user has no projects assigned, return an empty list instead of error
                    if not user_projects:
                        return jsonify([])

                    allowed_project_ids = {
                        pid for pid, perms in user_projects.items()
                        if perms["can_read"] or perms["can_write"] or perms["is_owner"]
                    }

                    filtered_projects = [
                        {
                            "id": p["id"],
                            "name": p["name"],
                            "lastModified": p["lastModified"],
                            "thumbnail": p["thumbnail"],
                            "permissions": user_projects[p["id"]],
                            "owner": get_project_owners(p["id"]),
                            "readme": p["readme"],
                        }
                        for p in active_projects
                        if p["id"] in allowed_project_ids
                    ]
                else:
                    # No auth, return all
                    filtered_projects = [
                        {
                            "id": p["id"],
                            "name": p["name"],
                            "lastModified": p["lastModified"],
                            "thumbnail": p["thumbnail"],
                            "permissions": None,
                            "owner": [],
                            "readme": p["readme"],
                        }
                        for p in active_projects
                    ]

                return jsonify(filtered_projects)

            except Exception as e:
                logger.exception(f"Error retrieving projects: {e}")
                return jsonify({"error": str(e)}), 500

        logger.info("Route registered: /projects")

    except Exception as e:
        logger.exception(f"Error registering routes: {e}")
        raise  # Re-raise to be handled by the parent function

