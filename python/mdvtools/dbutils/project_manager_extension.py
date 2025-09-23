from flask import Flask, jsonify, request
import os
import shutil
import zipfile
import io
from mdvtools.server_extension import MDVServerExtension
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprintProtocol
from mdvtools.dbutils.dbservice import ProjectService, UserProjectService
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


class ProjectManagerExtension(MDVServerExtension):
    def register_global_routes(self, app: Flask, config: dict):
        ENABLE_AUTH = app.config.get('ENABLE_AUTH', False)
        if ENABLE_AUTH:
            from flask import session
            from mdvtools.auth.authutils import user_cache

        @app.route("/create_project", methods=["POST"])
        def create_project():   
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
        def import_project():
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

                    # Reject entries with absolute paths or “..”
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
                p.set_editable(True)
                p.serve(app=app, open_browser=False, backend_db=True)
                

                
                # Initialize the project and register it using project name if valid
                if project_name is not None:
                    new_project = ProjectService.add_new_project(path=project_path, name=project_name)
                else:
                    new_project = ProjectService.add_new_project(path=project_path)

                if new_project:
                    if ENABLE_AUTH:
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

    def register_project_routes(self, project: MDVProject, project_bp: ProjectBlueprintProtocol):
        pass

    def mutate_state_json(self, state_json: dict, project: MDVProject, app: Flask):
        pass

    def get_session_config(self):
        return {
            "createProject": False,
            "importProject": False
        }