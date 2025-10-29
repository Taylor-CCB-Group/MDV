from flask import Flask, jsonify, request, send_file, Response
from typing import Union, Tuple
import os
import shutil
import zipfile
import io
import tempfile
from mdvtools.server_extension import MDVProjectServerExtension
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprintProtocol
from mdvtools.dbutils.dbservice import ProjectService, UserProjectService
from mdvtools.dbutils.dbmodels import User
from mdvtools.auth.authutils import update_cache
from mdvtools.dbutils.routes import get_project_thumbnail
from mdvtools.logging_config import get_logger

from mdvtools.project_router import ProjectBlueprint

logger = get_logger(__name__)


REQUIRED_FILES = {"views.json", "state.json", "datasources.json"}

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


class ProjectManagerExtension(MDVProjectServerExtension):
    def register_global_routes(self, app: Flask, config: dict):
        ENABLE_AUTH = app.config.get('ENABLE_AUTH', False)
        session = None
        user_cache = None
        active_projects_cache = []
        user_project_cache = {}
        all_users_cache = []
        
        if ENABLE_AUTH:
            from flask import session
            from mdvtools.auth.authutils import user_cache, active_projects_cache, user_project_cache, all_users_cache

        @app.route("/create_project", methods=["POST"])
        def create_project() -> Union[Response, Tuple[Response, int]]:   
            """
            Creates a new project and updates the caches and database accordingly.
            """
            project_path = None
            next_id = None
            try:
                
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
                    
                # Create a new Project record in the database with the path
                logger.info("Adding new project to the database")
                new_project = ProjectService.add_new_project(path=project_path)

                if new_project:
                    # Step 5: Associate the admin user with the new project and grant all permissions
                    if ENABLE_AUTH:
                        assert session is not None
                        assert user_cache is not None
                        user_data = session.get('user')
                        if not user_data:
                            raise ValueError("User not found in session.")
                        current_user_id = user_data['id']
                        UserProjectService.add_or_update_user_project(
                            user_id=current_user_id,
                            project_id=new_project.id,
                            is_owner=True
                        )
                        auth_id = user_data["auth_id"]
                        owner_email = user_cache.get(auth_id, {}).get("email", "unknown")
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
                                "thumbnail": thumbnail,
                                "owner": [owner_email]
                            }
                        )
                        
                    # Return the new project info
                    return jsonify({
                        "id": new_project.id,
                        "name": new_project.name,
                        "status": "success"
                    })
                else:
                    logger.error("Failed to create project in database")
                    return jsonify({"error": "Failed to create project in database"}), 500

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

        @app.route("/import_project", methods=["POST"])
        def import_project() -> Union[Response, Tuple[Response, int]]:
            try:
                # Check if the request contains a file
                if 'file' not in request.files:
                    logger.error("In register_routes /import_project: Error - No project archive provided")
                    return jsonify({"error": "No project archive provided"}), 400
                
                project_file = request.files['file']
                project_name = request.form.get('name')

                # Get next available project ID
                next_id = ProjectService.get_next_project_id()
                project_path = os.path.join(app.config['projects_base_dir'], str(next_id))
                os.makedirs(project_path, exist_ok=True)
                
                file_stream = io.BytesIO(project_file.read())

                with zipfile.ZipFile(file_stream) as zf:
                    names = zf.namelist()

                    # Reject entries with absolute paths or ".."
                    bad = [n for n in names if n.startswith(("/", "\\")) or ".." in n]
                    if bad:
                        logger.error("In register_routes /import_project: Error - Unsafe entries in ZIP")
                        return jsonify({"error": "Invalid ZIP file: unsafe paths detected"}), 400

                    # Find the root directory of the mdv project
                    # TODO: when we cut & pasted this block from routes.py, we didn't copy the find_root_prefix function
                    # so that needs to be moved somewhere we can use it here
                    root = find_root_prefix(names)
                    if root is None:
                        logger.error("In register_routes /import_project: Error - Not a valid MDV project")
                        return jsonify({"error": "Not a valid MDV project"}), 400

                    # Select the files based in the mdv project
                    members = [n for n in names if n.startswith(root)]

                    # Extract them to the project path
                    zf.extractall(path=project_path, members=members)

                # If everything was under a single top-level folder, flatten it
                if root:
                    subdir = os.path.join(project_path, root.rstrip("/"))
                    for item in os.listdir(subdir):
                        shutil.move(os.path.join(subdir, item), project_path)
                    os.rmdir(subdir)
                
                # # Create a new MDV project out of the new path and files copied
                p = MDVProject(project_path, backend_db=True)
                # Respect permission in imported state.json, default to non-editable if unspecified
                try:
                    state = p.state or {}
                    perm = (state.get('permission') or '').lower()
                    is_editable = True if perm == 'edit' else False if perm == 'view' else False
                    p.set_editable(is_editable)
                except Exception:
                    p.set_editable(False)
                p.serve(app=app, open_browser=False, backend_db=True)
                

                
                # Initialize the project and register it using project name if valid
                if project_name is not None:
                    new_project = ProjectService.add_new_project(path=project_path, name=project_name)
                else:
                    new_project = ProjectService.add_new_project(path=project_path)

                if new_project:
                    if ENABLE_AUTH:
                        assert session is not None
                        assert user_cache is not None
                        user_data = session.get("user")
                        if not user_data:
                            raise ValueError("User not found in session.")
                        current_user_id = user_data['id']
                        # Create a new entry with the new project and the current user
                        UserProjectService.add_or_update_user_project(
                            user_id=current_user_id,
                            project_id=new_project.id,
                            is_owner=True
                        )
                        auth_id = user_data["auth_id"]
                        owner_email = user_cache.get(auth_id, {}).get("email", "unknown")
                        thumbnail = get_project_thumbnail(project_path)

                        # Update cache with the new project data and user id
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
                                "thumbnail": thumbnail,
                                "owners": [owner_email]
                            }
                        )
                # Return the new project id and name
                logger.info("Import successfull. Returning the success response.")
                return jsonify({
                    "id": new_project.id,
                    "name": new_project.name,
                    "status": "success"
                })
                
            except Exception as e:
                logger.exception(f"In register_routes - /import_project : Error importing project: {e}")
                logger.info("started rollabck")
                # Clean up on error
                project_path = locals().get('project_path')
                next_id = locals().get('next_id')
                
                # Clean up project directory if it was created
                if project_path and os.path.exists(project_path):
                    try:
                        shutil.rmtree(project_path)
                        logger.info("In register_routes -/import_project : Rolled back project directory creation as db entry is not added")
                    except Exception as cleanup_error:
                        logger.exception(f"In register_routes -/import_project : Error during cleanup: {cleanup_error}")
                
                # Remove from blueprints if registered
                if next_id is not None and str(next_id) in ProjectBlueprint.blueprints:
                    try:
                        del ProjectBlueprint.blueprints[str(next_id)]
                        logger.info("In register_routes -/import_project : Rolled back ProjectBlueprint.blueprints as db entry is not added")
                    except Exception as blueprint_error:
                        logger.exception(f"In register_routes -/import_project : Error removing blueprint: {blueprint_error}")
                return jsonify({"error": str(e)}), 500
            
        logger.info("Route registered: /import_project")

        @app.route("/export_project/<int:project_id>", methods=["GET"])
        def export_project(project_id: int) -> Union[Response, Tuple[Response, int]]:
            try:

                # Fetch the project using the provided project id
                project = ProjectService.get_project_by_id(project_id)

                # No project with the given project id found
                if project is None:
                    logger.error(f"In register_routes - /export_project Error: Project with ID {project_id} not found in database")
                    return jsonify({"error": f"Project with ID {project_id} not found in database"})

                if ENABLE_AUTH:
                    if session is None:
                        raise ValueError("Session not available.")
                    user = session.get('user')
                    if not user:
                        raise ValueError("User not found in session.")
                    user_id = user["id"]
                    user_projects = user_project_cache.get(user_id) if user_project_cache is not None else None
                    if not user_projects or not user_projects.get(int(project_id), {}).get("is_owner", False):
                        logger.error(f"User does not have ownership of project {project_id}")
                        return jsonify({"error": "Only the project owner can export the project."}), 403
                    
                if project.access_level != 'editable':
                    logger.error(f"Project with ID {project_id} is not editable.")
                    return jsonify({"error": "This project is not editable and cannot be exported."}), 403
                    
                if project.path is None:
                    logger.error(f"In register_routes - /export_project Error: Project with ID {project_id} has no path in database")
                    return jsonify({"error": f"Project with ID {project_id} has no path in the database"})
                
                if project.name is None:
                    project_name = "unnamed_project"
                else:
                    project_name = project.name
                                
                # Create a temporary directory
                with tempfile.TemporaryDirectory() as temp_dir:
                    file_name = f"{project_name}"
                    file_path = os.path.join(temp_dir, file_name)
                    
                    # Create an archive from the project path
                    logger.info(f"Creating archive for project {project_name} at path: {project.path}...")
                    zip_path = shutil.make_archive(
                        file_path,
                        "zip",
                        project.path
                    )
                    logger.info(f"Archive created at: {zip_path}")
                    # Return the zip file
                    return send_file(
                        path_or_file=zip_path,
                        mimetype="application/zip",
                        as_attachment=True,
                        download_name=f"{project_name}.mdv.zip"
                    )

            except Exception as e:
                logger.exception(f"In register_routes - /export_project : Unexpected error while exporting project with project id - '{project_id}': {e}")
                return jsonify({"error": str(e)}), 500
        
        print("Route registered: /export_project/<project_id>")

        @app.route("/delete_project/<int:project_id>", methods=["DELETE"])
        def delete_project(project_id: int) -> Union[Response, Tuple[Response, int]]:
            #project_removed_from_blueprints = False
            nonlocal active_projects_cache
            try:
                logger.info(f"Deleting project '{project_id}'")
                
                # Step 2: Check if the user is owner using in-memory cache
                if ENABLE_AUTH:
                    if session is None:
                        raise ValueError("Session not available.")
                    user = session.get('user')
                    if not user:
                        raise ValueError("User not found in session.")
                    user_id = user["id"]
                    user_projects = user_project_cache.get(user_id) if user_project_cache is not None else None
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
        def rename_project(project_id: int) -> Union[Response, Tuple[Response, int]]:
            # Retrieve the new project name from the multipart/form-data payload
            new_name = request.form.get("name")
            
            if not new_name:
                return jsonify({"error": "New name not provided"}), 400
            
            try:
                # Step 1: Check if project exists in active cache
                project = ProjectService.get_project_by_id(project_id)
                if not project:
                    logger.error(f"Project with ID {project_id} not found in db")
                    return jsonify({"error": f"Project with ID {project_id} not found"}), 404

                # Step 2: Check ownership using cache
                if ENABLE_AUTH:
                    if session is None:
                        raise ValueError("Session not available.")
                    user = session.get('user')
                    if not user:
                        raise ValueError("User not found in session.")
                    user_id = user["id"]
                    user_projects = user_project_cache.get(user_id) if user_project_cache is not None else None
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
        def change_project_access(project_id: int) -> Union[Response, Tuple[Response, int]]:
            """API endpoint to change the access level of a project (editable or read-only)."""
            try:
                # Get the new access level from the request
                new_access_level = request.form.get("type")

                # Validate the new access level
                if new_access_level not in ["read-only", "editable"]:
                    return jsonify({"error": "Invalid access level. Must be 'read-only' or 'editable'."}), 400
            

                # Step 3: Check ownership from the user_project_cache
                if ENABLE_AUTH:
                    if session is None:
                        raise ValueError("Session not available.")
                    user = session.get('user')
                    if not user:
                        raise ValueError("User not found in session.")
                    user_id = user["id"]
                    user_projects = user_project_cache.get(user_id) if user_project_cache is not None else None
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
        def share_project(project_id: int) -> Union[Response, Tuple[Response, int]]:
            """Fetch users with whom the project is shared along with their permissions."""
            try:
                logger.info(f"Sharing project '{project_id}'")

                # Step 1: Skip authentication if not enabled
                if not ENABLE_AUTH:
                    logger.info("Authentication is disabled, skipping authentication check.")
                    return jsonify({"error": "Authentication is disabled, no action taken."})
                
                if session is None:
                    raise ValueError("Session not available.")
                user = session.get('user')
                if not user:
                    raise ValueError("User not found in session.")
                user_id = user["id"]

                user_permissions = user_project_cache.get(user_id, {}).get(int(project_id)) if user_project_cache is not None else {}
                if not user_permissions or not user_permissions.get("is_owner"):
                    return jsonify({"error": "Only the project owner can share the project"}), 403
                

                # Step 3: Compile shared users list
                shared_users_list = []
                for uid, permissions in user_project_cache.items():
                    proj_perm = permissions.get(project_id)
                    if proj_perm:
                        user_data = user_cache.get(uid) if user_cache is not None else None
                        if not user_data:
                            # Fallback: Try fetching from the DB if not in cache
                            user_obj = User.query.get(uid)
                            if not user_obj:
                                continue  # Skip if user doesn't exist in DB
                            user_data = {
                                "id": user_obj.id,
                                "email": user_obj.email,
                                "auth_id": user_obj.auth_id,
                                "is_admin": user_obj.is_admin
                            }
                            # Optionally update cache to avoid this hit next time
                            if user_cache is not None:
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
        def add_user_to_project(project_id: int) -> Union[Response, Tuple[Response, int]]:
            """Add a user to the project with specified permissions."""
            try:
                logger.info(f"Adding user to project '{project_id}'")

                # Step 1: Skip authentication if not enabled
                if not ENABLE_AUTH:
                    logger.info("Authentication is disabled, skipping authentication check.")
                    return jsonify({"error": "Authentication is disabled, no action taken."})

                if session is None:
                    raise ValueError("Session not available.")
                user = session.get('user')
                if not user:
                    raise ValueError("User not found in session.")
                user_id = user["id"]

                # Step 2: Check if current user is owner of the project
                user_permissions = user_project_cache.get(user_id, {}).get(int(project_id)) if user_project_cache is not None else {}
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
        def edit_user_permission(project_id: int, user_id: int) -> Union[Response, Tuple[Response, int]]:
            """Edit user permissions for a project."""
            try:
                logger.info(f"Editing permissions for user '{user_id}' in project '{project_id}'")

                # Step 1: Ensure authentication is enabled
                if not ENABLE_AUTH:
                    logger.info("Authentication is disabled, skipping authentication check.")
                    # If authentication is disabled, simply return and stop execution
                    return jsonify({"Error": "Authentication is disabled, no action taken."})

                if session is None:
                    raise ValueError("Session not available.")
                user = session.get('user')
                if not user:
                    raise ValueError("User not found in session.")
                current_user_id = user["id"]
                
                # Step 2: Validate if current user is the owner
                user_permissions = user_project_cache.get(current_user_id, {}).get(int(project_id)) if user_project_cache is not None else {}
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
                if user_project_cache is not None:
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
        def delete_user_from_project(project_id: int, user_id: int) -> Union[Response, Tuple[Response, int]]:
            """Remove a user from the project."""
            try:
                logger.info(f"Removing user '{user_id}' from project '{project_id}'")

                # Step 1: Ensure authentication is enabled
                if not ENABLE_AUTH:
                    logger.info("Authentication is disabled, skipping authentication check.")
                    # If authentication is disabled, simply return and stop execution
                    return jsonify({"error": "Authentication is disabled, no action taken."})

                if session is None:
                    raise ValueError("Session not available.")
                user = session.get('user')
                if not user:
                    raise ValueError("User not found in session.")
                current_user_id = user["id"]

                # Step 2: Validate ownership
                user_permissions = user_project_cache.get(current_user_id, {}).get(int(project_id)) if user_project_cache is not None else {}
                if not user_permissions or not user_permissions.get("is_owner"):
                    return jsonify({"error": "Only the project owner can remove users"}), 403

                # Step 4: Remove the user from the project using the service method
                UserProjectService.remove_user_from_project(user_id, project_id)

                # Step 5: Check if the user is part of the project
                # Step 5: Check if the user is part of the project
                user_permissions_cache = user_project_cache.get(user_id) if user_project_cache is not None else None
                if user_permissions_cache:
                    if user_project_cache is not None and user_id in user_project_cache and project_id in user_project_cache[user_id]:
                        del user_project_cache[user_id][project_id]
                    else:
                        logger.warning(f"Project {project_id} not found in cache for user {user_id}")
                return jsonify({"message": "User removed successfully"}), 200
            
            except Exception as e:
                logger.exception(f"Error in delete_user_from_project: {e}")
                return jsonify({"error": str(e)}), 500
        logger.info("Route registered: /projects/<int:project_id>/share/<int:user_id>/delete")
        
    def register_routes(self, project: MDVProject, project_bp: ProjectBlueprintProtocol):
        pass

    def mutate_state_json(self, state_json: dict, project: MDVProject, app: Flask):
        pass

    def get_session_config(self):
        return {
            "createProject": True,
            "importProject": True,
            "exportProject": True,
            "deleteProject": True,
            "renameProject": True,
            "changeProjectAccess": True,
            "shareProject": True,
            "editUserPermissions": True,
            "removeUserFromProject": True
        }