import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def register_routes(app, ENABLE_AUTH):
    from flask import abort, request, jsonify, redirect, url_for, render_template
    from mdvtools.auth.authutils import maybe_require_user, update_cache, active_projects_cache, user_project_cache,user_cache, all_users_cache
    from mdvtools.dbutils.mdv_server_app import serve_projects_from_filesystem
    import os
    import shutil
    from mdvtools.mdvproject import MDVProject
    from mdvtools.project_router import ProjectBlueprint_v2 as ProjectBlueprint
    from mdvtools.dbutils.dbmodels import User
    from mdvtools.dbutils.dbservice import ProjectService, UserProjectService

    """Register routes with the Flask app."""
    logger.info("Registering routes...")
    # if len(active_projects_cache) == 0:
    #     logger.info("active_projects_cache is empty, calling cache_user_projects()")
    #     with app.app_context():
    #         cache_user_projects() # this won't work in localhost
    #         logger.info(f"active_projects_cache populated with {len(active_projects_cache)} projects")

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

        @app.route('/rescan_projects')
        @maybe_require_user(ENABLE_AUTH)
        def rescan_projects(user):
            if ENABLE_AUTH:
                is_admin = user.get("is_admin", False)
                if not is_admin:
                    abort(403)  # Forbidden

            #Serve the project
            serve_projects_from_filesystem(app, app.config["projects_base_dir"])
            return redirect(url_for('index'))


        @app.route('/projects')
        @maybe_require_user(ENABLE_AUTH)  # Pass ENABLE_AUTH here
        def get_projects(user):
            """
            Fetches the list of active projects from cache or database depending on the authentication.
            """
            logger.info(" /projects queried...")

            try:
                # Step 1: Get projects directly from cache
                if ENABLE_AUTH:
                    active_projects = active_projects_cache
                else:
                    active_projects = ProjectService.get_active_projects()
                # logger.info(f"Active projects from cache: {active_projects}")

                if ENABLE_AUTH:
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
                            "thumbnail": p["thumbnail"]
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
                            "thumbnail": p["thumbnail"]
                        }
                        for p in active_projects
                    ]

                return jsonify(filtered_projects)

            except Exception as e:
                logger.exception(f"Error retrieving projects: {e}")
                return jsonify({"error": str(e)}), 500

        logger.info("Route registered: /projects")


        @app.route("/create_project", methods=["POST"])
        @maybe_require_user(ENABLE_AUTH)  # Pass ENABLE_AUTH here to handle authentication
        def create_project(user):   
            """
            Creates a new project and updates the caches and database accordingly.
            """
            project_path = None
            next_id = None
            try:
                user_data = user  # Retrieve user data from the decorated function

                # Check if authentication is enabled and user is admin
                # Check if authentication is enabled and user is admin
                #if ENABLE_AUTH:
                #    # Ensure the user is an admin before allowing project creation
                #    if not user_data.get("is_admin", False):
                #        return jsonify({"error": "Only admins can create projects."}), 403

                logger.info("Creating project")
                
                # Get the next available ID
                next_id = ProjectService.get_next_project_id()
                if next_id is None:
                    logger.error("In register_routes: Error- Failed to determine next project ID from db")
                    return jsonify({"error": "Failed to determine next project ID from db"}), 500

                # Create the project directory path
                project_path = os.path.join(app.config['projects_base_dir'], str(next_id))

                # Create and serve the MDVProject
                try:
                    logger.info("Creating and serving the new project")
                    p = MDVProject(project_path, backend_db= True)
                    p.set_editable(True)
                    p.serve(app=app, open_browser=False, backend_db=True)
                except Exception as e:
                    logger.exception(f"In register_routes: Error serving MDVProject: {e}")
                    return jsonify({"error": "Failed to serve MDVProject"}), 500

                def get_project_thumbnail(project_path):
                    """Extract the first available viewImage from a project's views."""
                    try:
                        mdv_project = MDVProject(project_path)
                        return next((v["viewImage"] for v in mdv_project.views.values() if "viewImage" in v), None)
                    except Exception as e:
                        logger.exception(f"Error extracting thumbnail for project at {project_path}: {e}")
                        return None
                    
                # Create a new Project record in the database with the path
                logger.info("Adding new project to the database")
                new_project = ProjectService.add_new_project(path=project_path)

                if new_project:
                    # Step 5: Associate the admin user with the new project and grant all permissions
                    if ENABLE_AUTH:
                        current_user_id = user_data['id']
                        UserProjectService.add_or_update_user_project(
                            user_id=current_user_id,
                            project_id=new_project.id,
                            is_owner=True
                        )

                        # Generate thumbnail
                        thumbnail = get_project_thumbnail(project_path)
                        
                        # Step 6: Update caches for the admin user and new project
                        update_cache(
                            user_id=current_user_id,
                            project_id=new_project.id,
                            permissions={
                                "can_read": True,
                                "can_write": True,
                                "is_owner": True
                            },
                            project_data={
                                "id": new_project.id,
                                "name": new_project.name,
                                "lastModified": new_project.update_timestamp.strftime("%Y-%m-%d %H:%M:%S"),
                                "thumbnail": thumbnail
                            }
                        )
                        
                    # Return the new project info
                    return jsonify({
                        "id": new_project.id,
                        "name": new_project.name,
                        "status": "success"
                    })     

            except Exception as e:
                logger.error(f"In register_routes - /create_project : Error creating project: {e}")
                logger.info("started rollabck")
                # Rollback: Clean up the projects filesystem directory if it was created
                if project_path and os.path.exists(project_path):
                    try:
                        shutil.rmtree(project_path)
                        logger.info("In register_routes -/create_project : Rolled back project directory creation as db entry is not added")
                    except Exception as cleanup_error:
                        logger.exception(f"In register_routes -/create_project : Error during rollback cleanup: {cleanup_error}")

                # Optional: Remove project routes from Flask app if needed
                if next_id is not None and str(next_id) in ProjectBlueprint.blueprints:
                    del ProjectBlueprint.blueprints[str(next_id)]
                    logger.info("In register_routes -/create_project : Rolled back ProjectBlueprint.blueprints as db entry is not added")
                
                return jsonify({"error": str(e)}), 500

        logger.info("Route registered: /create_project")

        @app.route("/delete_project/<int:project_id>", methods=["DELETE"])
        @maybe_require_user(ENABLE_AUTH)
        def delete_project(user, project_id: int):
            #project_removed_from_blueprints = False
            nonlocal active_projects_cache
            try:
                logger.info(f"Deleting project '{project_id}'")

                user_id = user["id"] if ENABLE_AUTH else None

    
                # Step 2: Check if the user is owner using in-memory cache
                if ENABLE_AUTH:
                    user_projects = user_project_cache.get(user_id)
                    if not user_projects or not user_projects.get(int(project_id), {}).get("is_owner", False):
                        logger.error(f"User does not have ownership of project {project_id}")
                        return jsonify({"error": "Only the project owner can delete the project."}), 403


                # Step 3: Get full project object for access_level check
                project_obj = ProjectService.get_project_by_id(project_id)
                if not project_obj:
                    logger.error(f"Project with ID {project_id} not found in database")
                    return jsonify({"error": f"Project with ID {project_id} not found in database"}), 404

                # Step 4: Check if project is editable
                if project_obj.access_level != 'editable':
                    logger.error(f"Project with ID {project_id} is not editable.")
                    return jsonify({"error": "This project is not editable and cannot be deleted."}), 403


                # Remove the project from the ProjectBlueprint.blueprints dictionary
                if str(project_id) in ProjectBlueprint.blueprints:
                    del ProjectBlueprint.blueprints[str(project_id)]
                    #project_removed_from_blueprints = True  # Mark as removed
                    logger.info(f"In register_routes - /delete_project : Removed project '{project_id}' from ProjectBlueprint.blueprints")
                
                # Soft delete the project
                delete_status = ProjectService.soft_delete_project(project_id)
                if not delete_status:
                    logger.error(f"In delete_project: Error - Failed to soft delete project {project_id} in the database")
                    return jsonify({"error": "Failed to soft delete project in db"}), 500

                # Step 7: Remove from cache
                if ENABLE_AUTH:
                    
                    active_projects_cache = [p for p in active_projects_cache if p["id"] != project_id]
                    logger.info(f"Removed project '{project_id}' from in-memory cache")

                # Step 8: Return success response
                return jsonify({"message": f"Project '{project_id}' deleted successfully."})

            except Exception as e:
                logger.exception(f"In register_routes - /delete_project: Error deleting project '{project_id}': {e}")
                return jsonify({"error": str(e)}), 500

        logger.info("Route registered: /delete_project/<project_id>")

        @app.route("/projects/<int:project_id>/rename", methods=["PUT"])
        @maybe_require_user(ENABLE_AUTH)
        def rename_project(user, project_id: int):
            # Retrieve the new project name from the multipart/form-data payload
            new_name = request.form.get("name")
            
            if not new_name:
                return jsonify({"error": "New name not provided"}), 400
            
            try:
                user_id = user["id"] if ENABLE_AUTH else None

                # Step 1: Check if project exists in active cache
                project = ProjectService.get_project_by_id(project_id)
                if not project:
                    logger.error(f"Project with ID {project_id} not found in db")
                    return jsonify({"error": f"Project with ID {project_id} not found"}), 404

                # Step 2: Check ownership using cache
                if ENABLE_AUTH:
                    user_projects = user_project_cache.get(user_id)
                    if not user_projects or not user_projects.get(int(project_id), {}).get("is_owner", False):
                        logger.error(f"User does not have ownership of project {project_id}")
                        return jsonify({"error": "Only the project owner can rename the project."}), 403

        
                # Step 4: Check if project is editable
                if project.access_level != 'editable':
                    logger.error(f"Project with ID {project_id} is not editable.")
                    return jsonify({"error": "This project is not editable and cannot be renamed."}), 403

                # Attempt to rename the project
                rename_status = ProjectService.update_project_name(project_id, new_name)

                if not rename_status:
                    logger.error(f"In register_routes - /rename_project Error: The project with ID '{project_id}' not found in db")
                    return jsonify({"error": f"Failed to rename project '{project_id}' in db"}), 500
                
                # Step 6: Update the name in active_projects_cache
                if ENABLE_AUTH:
                    for proj in active_projects_cache:
                        if proj["id"] == project_id:
                            proj["name"] = new_name
                            break

                logger.info(f"In rename_project: Updated cache for active projects, renamed project '{project_id}'.")

                # Step 7: Return success response
                return jsonify({"status": "success", "id": project_id, "new_name": new_name}), 200

            except Exception as e:
                logger.exception(f"In register_routes - /rename_project : Error renaming project '{project_id}': {e}")
                return jsonify({"error": str(e)}), 500

        logger.info("Route registered: /projects/<int:project_id>/rename")

        @app.route("/projects/<int:project_id>/access", methods=["PUT"])
        @maybe_require_user(ENABLE_AUTH)
        def change_project_access(user, project_id: int):
            """API endpoint to change the access level of a project (editable or read-only)."""
            try:
                # Get the new access level from the request
                new_access_level = request.form.get("type")

                # Validate the new access level
                if new_access_level not in ["read-only", "editable"]:
                    return jsonify({"error": "Invalid access level. Must be 'read-only' or 'editable'."}), 400
            
                user_id = user["id"] if ENABLE_AUTH else None

                # Step 3: Check ownership from the user_project_cache
                if ENABLE_AUTH:
                    user_projects = user_project_cache.get(user_id)
                    if not user_projects or not user_projects.get(int(project_id), {}).get("is_owner", False):
                        logger.error(f"User does not have ownership of project {project_id}")
                        return jsonify({"error": "Only the project owner can change the access level."}), 403


                # Call the service method to change the access level
                access_level, message, status_code = ProjectService.change_project_access(project_id, new_access_level)

                if access_level is None:
                    return jsonify({"error": message}), status_code

                return jsonify({"status": "success", "access_level": access_level}), 200

            except Exception as e:
                logger.exception(f"In register_routes - /access : Unexpected error while changing access level for project '{project_id}': {e}")
                return jsonify({"error": "An unexpected error occurred."}), 500
        
        logger.info("Route registered: /projects/<int:project_id>/access")

        @app.route("/projects/<int:project_id>/share", methods=["GET"])
        @maybe_require_user(ENABLE_AUTH)
        def share_project(user, project_id: int):
            """Fetch users with whom the project is shared along with their permissions."""
            try:
                logger.info(f"Sharing project '{project_id}'")

                # Step 1: Skip authentication if not enabled
                if not ENABLE_AUTH:
                    logger.info("Authentication is disabled, skipping authentication check.")
                    return jsonify({"error": "Authentication is disabled, no action taken."})
                
                user_id = user["id"] if ENABLE_AUTH else None

                user_permissions = user_project_cache.get(user_id, {}).get(int(project_id))
                if not user_permissions or not user_permissions.get("is_owner"):
                    return jsonify({"error": "Only the project owner can share the project"}), 403
                

                # Step 3: Compile shared users list
                shared_users_list = []
                for uid, permissions in user_project_cache.items():
                    proj_perm = permissions.get(project_id)
                    if proj_perm:
                        user_data = user_cache.get(uid)
                        if not user_data:
                            # Fallback: Try fetching from the DB if not in cache
                            user_obj = User.query.get(uid)
                            if not user_obj:
                                continue  # Skip if user doesn't exist in DB
                            user_data = {
                                "id": user_obj.id,
                                "email": user_obj.email,
                                "auth0_id": user_obj.auth0_id,
                                "is_admin": user_obj.is_admin
                            }
                            # Optionally update cache to avoid this hit next time
                            user_cache[uid] = user_data
                            if user_data not in all_users_cache:
                                all_users_cache.append(user_data)

                        shared_users_list.append({
                            "id": user_data["id"],
                            "email": user_data["email"],
                            "permission": (
                                "Owner" if proj_perm["is_owner"] else "Edit" if proj_perm["can_write"] else "View"
                            )
                        })


                # Step 4: Get unshared users for dropdown
                shared_user_ids = {u["id"] for u in shared_users_list}
                all_users = [
                    {"id": u["id"], "email": u["email"]}
                    for u in all_users_cache
                    if u["id"] not in shared_user_ids
                ]
                
                # Return the list of users with permissions, and all users for the dropdown
                return jsonify({
                    "shared_users": shared_users_list,
                    "all_users": all_users  # List of all users to populate the dropdown
                })
            
            except Exception as e:
                logger.exception(f"Error in share_project: {e}")
                return jsonify({"error": str(e)}), 500
        logger.info("Route registered: /projects/<int:project_id>/share- GET")
            
        @app.route("/projects/<int:project_id>/share", methods=["POST"])
        @maybe_require_user(ENABLE_AUTH)
        def add_user_to_project(user, project_id: int):
            """Add a user to the project with specified permissions."""
            try:
                logger.info(f"Adding user to project '{project_id}'")

                # Step 1: Skip authentication if not enabled
                if not ENABLE_AUTH:
                    logger.info("Authentication is disabled, skipping authentication check.")
                    return jsonify({"error": "Authentication is disabled, no action taken."})

                user_id = user["id"] if ENABLE_AUTH else None

                # Step 2: Check if current user is owner of the project
                user_permissions = user_project_cache.get(user_id, {}).get(int(project_id))
                if not user_permissions or not user_permissions.get("is_owner"):
                    return jsonify({"error": "Only the project owner can share the project"}), 403


                # Step 5: Get data from the POST request
                data = request.get_json()
                target_user_id = data.get('user_id')  # User to be added
                permission = data.get('permission', 'view')  # "view", "edit", or "owner"

                if not target_user_id or permission not in ["view", "edit", "owner"]:
                    return jsonify({"error": "Invalid user or permission"}), 400
                
                # Step 4: Determine permission flags
                is_owner = permission == "owner"
                can_write = permission in ["edit", "owner"]
                #can_read = True  # Always true if added to a project

                # Step 5: Update the DB via service
                UserProjectService.add_or_update_user_project(
                    user_id=target_user_id,
                    project_id=project_id,
                    is_owner=is_owner,
                    can_write=can_write
                )
                
                # Update user_project_cache
                if target_user_id not in user_project_cache:
                    user_project_cache[target_user_id] = {}

                user_project_cache[target_user_id][project_id] = {
                    "can_read": (permission in ["view", "edit", "owner"]),
                    "can_write": (permission in ["edit", "owner"]),
                    "is_owner": (permission == "owner")
                }

                logger.info(f"User {target_user_id} added to project {project_id} with '{permission}' permission.")

                return jsonify({"message": f"User {target_user_id} added to project {project_id} with {permission} permission."}), 200

            except Exception as e:
                logger.exception(f"Error in add_user_to_project: {e}")
                return jsonify({"error": str(e)}), 500
        logger.info("Route registered: /projects/<int:project_id>/share- POST")


        @app.route("/projects/<int:project_id>/share/<int:user_id>/edit", methods=["POST"])
        @maybe_require_user(ENABLE_AUTH)
        def edit_user_permission(user, project_id: int, user_id: int):
            """Edit user permissions for a project."""
            try:
                logger.info(f"Editing permissions for user '{user_id}' in project '{project_id}'")

                # Step 1: Ensure authentication is enabled
                if not ENABLE_AUTH:
                    logger.info("Authentication is disabled, skipping authentication check.")
                    # If authentication is disabled, simply return and stop execution
                    return jsonify({"Error": "Authentication is disabled, no action taken."})

                current_user_id = user["id"] if ENABLE_AUTH else None # The ID of the authenticated user

                # Step 2: Validate if current user is the owner
                user_permissions = user_project_cache.get(current_user_id, {}).get(int(project_id))
                if not user_permissions or not user_permissions.get("is_owner"):
                    return jsonify({"error": "Only the project owner can edit permissions"}), 403


                # Step 3: Extract new permission from request
                assert request.json, "Request JSON is empty"
                new_permission = request.json.get("permission", "").lower()
                if new_permission not in ["view", "edit", "owner"]:
                    return jsonify({"error": "Invalid permission value"}), 400
                
                is_owner = new_permission == "owner"
                can_write = new_permission in ["edit", "owner"]
                can_read = True  # Always true for project access

                # Step 4: Update DB via service
                UserProjectService.add_or_update_user_project(
                    user_id=user_id,
                    project_id=project_id,
                    is_owner=is_owner,
                    can_write=can_write
                )
                
                # Step 5: Update cache
                if user_id not in user_project_cache:
                    user_project_cache[user_id] = {}

                user_project_cache[user_id][project_id] = {
                    "can_read": can_read,
                    "can_write": can_write,
                    "is_owner": is_owner
                }

                logger.info(f"Updated permissions for user {user_id} in project {project_id}: {new_permission}")
                return jsonify({"message": "Permissions updated successfully"}), 200
            
            except Exception as e:
                logger.exception(f"Error in edit_user_permission: {e}")
                return jsonify({"error": str(e)}), 500
        logger.info("Route registered: /projects/<int:project_id>/share/<int:user_id>/edit")
            
        @app.route("/projects/<int:project_id>/share/<int:user_id>/delete", methods=["POST"])
        @maybe_require_user(ENABLE_AUTH)
        def delete_user_from_project(user, project_id: int, user_id: int):
            """Remove a user from the project."""
            try:
                logger.info(f"Removing user '{user_id}' from project '{project_id}'")

                # Step 1: Ensure authentication is enabled
                if not ENABLE_AUTH:
                    logger.info("Authentication is disabled, skipping authentication check.")
                    # If authentication is disabled, simply return and stop execution
                    return jsonify({"error": "Authentication is disabled, no action taken."})

                current_user_id = user["id"] if ENABLE_AUTH else None # The ID of the authenticated user

                # Step 2: Validate ownership
                user_permissions = user_project_cache.get(current_user_id, {}).get(int(project_id))
                if not user_permissions or not user_permissions.get("is_owner"):
                    return jsonify({"error": "Only the project owner can remove users"}), 403

                # Step 4: Remove the user from the project using the service method
                UserProjectService.remove_user_from_project(user_id, project_id)

                # Step 5: Check if the user is part of the project
                # Step 5: Check if the user is part of the project
                user_permissions_cache = user_project_cache.get(user_id)
                if user_permissions_cache:
                    if project_id in user_project_cache[user_id]:
                        del user_project_cache[user_id][project_id]
                    else:
                        logger.warning(f"Project {project_id} not found in cache for user {user_id}")
                return jsonify({"message": "User removed successfully"}), 200
            
            except Exception as e:
                logger.exception(f"Error in delete_user_from_project: {e}")
                return jsonify({"error": str(e)}), 500
        logger.info("Route registered: /projects/<int:project_id>/share/<int:user_id>/delete")


    except Exception as e:
        logger.exception(f"Error registering routes: {e}")
        raise  # Re-raise to be handled by the parent function

