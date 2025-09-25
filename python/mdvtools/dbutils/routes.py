from mdvtools.logging_config import get_logger
logger = get_logger(__name__)

def register_routes(app, ENABLE_AUTH):
    from flask import abort, request, jsonify, session, redirect, url_for, render_template
    from mdvtools.auth.authutils import active_projects_cache, user_project_cache, user_cache, all_users_cache
    from mdvtools.dbutils.mdv_server_app import serve_projects_from_filesystem
    from mdvtools.dbutils.dbservice import ProjectService
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

        @app.route('/rescan_projects')
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
                serve_projects_from_filesystem(app, app.config["projects_base_dir"])
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

