import os
import time
import json
from sqlalchemy import text
from sqlalchemy.exc import OperationalError
from flask import Flask, render_template, jsonify
#from flask_sqlalchemy import SQLAlchemy
# import threading
# import random
# import string
# from flask import Flask, render_template, jsonify, request
from mdvtools.server import add_safe_headers
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprint
from mdvtools.dbutils.dbmodels import db, Project
#from mdvtools.dbutils.routes import register_global_routes
from mdvtools.dbutils.dbservice import ProjectService, FileService


def create_flask_app(config_name=None):
    """*Create and configure the Flask app."""
    app = Flask(__name__, template_folder='../templates', static_folder='/app/dist/flask')
    # this was causing a number of issues with various routes, changing this here seems to be the best way to fix it
    # as there isn't a clear single point of front-end that would consistently guarantee fixing it
    app.url_map.strict_slashes = False
    app.after_request(add_safe_headers)

    
    try:
        print("*Adding config.json details to app config")
        load_config(app, config_name)
    except Exception as e:
        print(f"Error loading configuration: {e}")
        exit(1)
    
    try:
        print("Creating base directory")
        create_base_directory(app)
    except Exception as e:
        print(f"Error creating base directory: {e}")
        exit(1)

    try:
        print("Initializing app with db")
        db.init_app(app)
    except Exception as e:
        print(f"Error initializing database: {e}")
        exit(1)

    try:
        print("Creating tables")
        with app.app_context():
            print("********* Waiting for DB to set up")
            wait_for_database()
            if not tables_exist():
                print("Creating database tables")
                db.create_all()
                print("************* Created database tables")
            else:
                print("Database tables already exist")

            # Routes registration and application setup
            print("Registering the blueprint (register_app)")
            ProjectBlueprint.register_app(app)
            

    except OperationalError as oe:
        print(f"OperationalError: {oe}")
        exit(1)
    except Exception as e:
        print(f"Error during app setup: {e}")
        exit(1)

    try:
        # Register routes
        print("Registering base routes: /, /projects, /create_project, /delete_project")
        register_routes(app)
    except Exception as e:
        print(f"Error registering routes: {e}")
        exit(1)

    return app



def wait_for_database():
    """Wait for the database to be ready before proceeding."""
    max_retries = 30
    delay = 5  # seconds

    for attempt in range(max_retries):
        try:
            # Test database connection using engine.connect()
            with db.engine.connect() as connection:
                connection.execute(text('SELECT 1'))
            print("*************** Database is ready! *************")
            return
        except OperationalError as oe:
            print(f"Database not ready, retrying in {delay} seconds... (Attempt {attempt + 1} of {max_retries})")
            time.sleep(delay)
        except Exception as e:
            print(f"An unexpected error occurred while waiting for the database: {e}")
            raise  # Re-raise the exception to be handled by the parent

    # If the loop completes without a successful connection
    error_message = "Error: Database did not become available in time."
    print(error_message)
    raise TimeoutError(error_message)



def load_config(app,config_name=None):
    try:
        config_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.json')
        with open(config_file_path) as config_file:
            config = json.load(config_file)
            #app.config['PREFERRED_URL_SCHEME'] = 'https'
            app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = config.get('track_modifications', False)
            app.config['upload_folder'] = config.get('upload_folder', '')
            app.config['projects_base_dir'] = config.get('projects_base_dir', '')
            app.config['db_host'] = config.get('db_container', '')
            print("Configuration loaded successfully!")
    except FileNotFoundError:
        print("Error: Configuration file not found.")
        raise
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        raise

    # Handle different environments
    try:
        if config_name == 'test':
            app.config['PREFERRED_URL_SCHEME'] = 'http'
            app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///:memory:'
        else:
            app.config['PREFERRED_URL_SCHEME'] = 'http'
            
            # Load sensitive data from Docker secrets
            def read_secret(secret_name):
                secret_path = f'/run/secrets/{secret_name}'
                try:
                    with open(secret_path, 'r') as secret_file:
                        return secret_file.read().strip()
                except FileNotFoundError as fnf_error:
                    print(f"Error: Secret '{secret_name}' not found. {fnf_error}")
                    raise  # Re-raise the exception to be handled by the parent

            db_user = read_secret('db_user')
            db_password = read_secret('db_password')
            db_name = read_secret('db_name')
            db_host = app.config['db_host']

            if not all([db_user, db_password, db_name, db_host]):
                raise ValueError("Error: One or more required secrets or configurations are missing.")
            
            app.config['SQLALCHEMY_DATABASE_URI'] = f'postgresql://{db_user}:{db_password}@{db_host}/{db_name}'

    except Exception as e:
        print(f"An unexpected error occurred while configuring the database: {e}")
        raise  # Re-raise the exception to be handled by the parent
       


# Function to create base directory if it doesn't exist
def create_base_directory(app):
    try:
        base_dir = app.config.get('projects_base_dir', 'mdv')
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
            print(f"Created base directory: {base_dir}")
        else:
            print(f"Base directory already exists: {base_dir}")
    except Exception as e:
        print(f'Function create_base_directory Error: {e}')
        raise

def tables_exist():
        inspector = db.inspect(db.engine)
        print("printing table names")
        print(inspector.get_table_names())
        return inspector.get_table_names()

def serve_projects_from_db(app):
    try:
        # Get all projects from the database
        print("Serving the projects present in both database and filesystem. Displaying the error if the path doesn't exist for a project")
        projects = Project.query.all()

        for project in projects:
            
            if project.is_deleted:
                print(f"Project with ID {project.id} is soft deleted.")
                continue

            if os.path.exists(project.path):
                try:
                    p = MDVProject(dir=project.path, id=str(project.id))
                    p.set_editable(True)
                    # todo: look up how **kwargs works and maybe have a shared app config we can pass around
                    p.serve(app=app, open_browser=False, backend=True)
                    print(f"Serving project: {project.path}")

                    # Update or add files in the database to reflect the actual files in the filesystem
                    for root, dirs, files in os.walk(project.path):
                        for file_name in files:
                            full_file_path = os.path.join(root, file_name)

                            # Use the utility function to add or update the file in the database
                            updated_file = FileService.add_or_update_file_in_project(
                                file_name=file_name,
                                file_path=full_file_path,
                                project_id=project.id
                            )
                            if updated_file is None:
                                print(f"Failed to add or update file '{file_name}' in the database.")
                            else:
                                print(f"Processed file in DB: {updated_file}")

                except Exception as e:
                    print(f"Error serving project '{project.path}': {e}")
                    raise
            else:
                print(f"Error: Project path '{project.path}' does not exist.")
                
    except Exception as e:
        print(f"Error serving projects from database: {e}")
        raise

def serve_projects_from_filesystem(app, base_dir):
    try:
        print("Serving the projects present in filesystem but missing in database")
        print(f"Scanning base directory: {base_dir}")

        # Get all project paths from the database
        projects_in_db = {project.path for project in Project.query.with_entities(Project.path).all()}
        print(f"Project paths in DB: {projects_in_db}")

        # Get all project directories in the filesystem
        project_paths_in_fs = {os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))}
        print(f"Project paths in filesystem: {project_paths_in_fs}")

        # Determine which project paths are in the filesystem but not in the database
        missing_project_paths = project_paths_in_fs - projects_in_db
        print(f"Missing project paths: {missing_project_paths}")

        # Iterate over missing project paths to create and serve them
        for project_path in missing_project_paths:
            print(f"Processing project path: {project_path}")
            
            if os.path.exists(project_path):
                try:
                    project_name = os.path.basename(project_path)

                    # Get the next ID from the database
                    next_id = db.session.query(db.func.max(Project.id)).scalar()
                    if next_id is None:
                        next_id = 1
                    else:
                        next_id += 1

                    p = MDVProject(dir=project_path,id= str(next_id))
                    p.set_editable(True)
                    p.serve(app=app, open_browser=False, backend=True) 
                    print(f"Serving project: {project_path}")

                    # Create a new Project record in the database with the default name
                    new_project = ProjectService.add_new_project(name=project_name, path=project_path)
                    if new_project is None:
                        raise ValueError(f"Failed to add project '{project_name}' to the database.")
                    else:
                        print(f"Added project to DB: {new_project}")
                    
                    # Add files from the project directory to the database
                    for root, dirs, files in os.walk(project_path):
                        for file_name in files:
                            # Construct the full file path
                            full_file_path = os.path.join(root, file_name)
                            
                            # Use the full file path when adding or updating the file in the database
                            new_file = FileService.add_or_update_file_in_project(
                                file_name=file_name,
                                file_path=full_file_path,
                                project_id=new_project.id
                            )
                            
                            if new_file is None:
                                raise ValueError(f"Failed to add file '{file_name}' to the database.")
                            else:
                                print(f"Added file to DB: {new_file}")


                except Exception as e:
                    print(f"In create_projects_from_filesystem: Error creating project at path '{project_path}': {e}")
                    raise
            else:
                print(f"In create_projects_from_filesystem: Error - Project path '{project_path}' does not exist.")
                
    except Exception as e:
        print(f"In create_projects_from_filesystem: Error retrieving projects from database: {e}")
        raise


def register_routes(app):
    """Register routes with the Flask app."""
    print("Registering routes...")

    try:
        @app.route('/')
        def index():
            try:
                return render_template('index.html')
            except Exception as e:
                print(f"Error rendering index: {e}")
                return jsonify({"status": "error", "message": str(e)}), 500

        print("Route registered: /")

        @app.route('/projects')
        def get_projects():
            print('/projects queried...')
            try:
                # Query the database to get all projects that aren't deleted
                projects = ProjectService.get_active_projects()

                # Return the list of projects with their IDs and names
                return jsonify([{"id": p.id, "name": p.name} for p in projects])
            except Exception as e:
                print(f"Error retrieving projects: {e}")
                return jsonify({"status": "error", "message": str(e)}), 500

        print("Route registered: /projects")

        @app.route("/create_project", methods=["POST"])
        def create_project():
            try:
                print("Creating project")
                
                # Get the next available ID
                next_id = ProjectService.get_next_project_id()
                if next_id is None:
                    print("Error: Failed to determine next project ID from db")
                    return jsonify({"status": "error", "message": "Failed to determine next project ID from db"}), 500

                # Create the project directory path
                project_path = os.path.join(app.config['projects_base_dir'], str(next_id))

                # Create and serve the MDVProject
                print("Creating and serving the new project")
                p = MDVProject(project_path)
                p.set_editable(True)
                p.serve(app=app, open_browser=False, backend=True)

                # Create a new Project record in the database with the path
                print("Adding new project to the database")
                new_project = ProjectService.add_new_project(path=project_path)

                if new_project:
                    return jsonify({"id": new_project.id, "name": new_project.name, "status": "success"})
                else:
                    print("Error: Failed to add new project to db")
                    return jsonify({"status": "error", "message": "Failed to add new project to db"}), 500

            except Exception as e:
                print(f"Error creating project: {e}")
                return jsonify({"status": "error", "message": str(e)}), 500

        print("Route registered: /create_project")

        @app.route("/delete_project/<project_id>", methods=["DELETE"])
        def delete_project(project_id: int):
            try:
                print(f"Deleting project '{project_id}'")
                
                # Find the project by ID 
                project = ProjectService.get_project_by_id(project_id)

                if project is None:
                    print(f"Error: Project with ID {project_id} not found in database")
                    return jsonify({"status": "error", "message": f"Project with ID {project_id} not found in database"}), 404

                # Remove the project from the ProjectBlueprint.blueprints dictionary
                if str(project_id) in ProjectBlueprint.blueprints:
                    del ProjectBlueprint.blueprints[str(project_id)]
                    print(f"Removed project '{project_id}' from ProjectBlueprint.blueprints")
                
                # Soft delete the project
                delete_status = ProjectService.soft_delete_project(project_id)

                if delete_status:
                    return jsonify({"status": "success"})
                else:
                    print("Error: Failed to soft delete project in db")
                    return jsonify({"status": "error", "message": "Failed to soft delete project in db"}), 500

            except Exception as e:
                print(f"Error deleting project '{project_id}': {e}")
                return jsonify({"status": "error", "message": str(e)}), 500

        print("Route registered: /delete_project/<project_id>")

    except Exception as e:
        print(f"Error registering routes: {e}")
        raise  # Re-raise to be handled by the parent function


"""
    
    print("Initialized app context")

        
    #db.create_all()
    #print("Created the database tables")

    #print("Registering the global routes")
    #register_global_routes(app, db, app.config['projects_base_dir'])
        
    print("Registering the blueprint(register_app)")
    ProjectBlueprint.register_app(app)
        
    print("Serve projects from database")
    serve_projects_from_db()

    print("Start- create_projects_from_filesystem")
    serve_projects_from_filesystem(app.config['projects_base_dir']) """

# Create the app object at the module level
app = create_flask_app()

with app.app_context():
    print("Serving projects from database")
    serve_projects_from_db(app)
    print("Starting - create_projects_from_filesystem")
    serve_projects_from_filesystem(app, app.config['projects_base_dir'])

if __name__ == '__main__':
    print("In main..")
    #wait_for_database()
    
    app.run(host='0.0.0.0', debug=True, port=5055)
    
    