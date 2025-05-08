import os
import time
import json
import logging
from sqlalchemy import text
from sqlalchemy.exc import OperationalError
from flask import Flask
from mdvtools.server import add_safe_headers
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprint_v2 as ProjectBlueprint
from mdvtools.dbutils.dbmodels import db, Project
from mdvtools.dbutils.routes import register_routes
from mdvtools.auth.register_auth_routes import register_auth_routes
from mdvtools.auth.authutils import register_before_request_auth, sync_auth0_users_to_db, cache_user_projects
from mdvtools.dbutils.dbservice import ProjectService, FileService

# Setup logging
logger = logging.getLogger(__name__)

#Read environment flag for authentication
ENABLE_AUTH = os.getenv("ENABLE_AUTH", "0").lower() in ["1", "true", "yes"]
logger.info(f"Authentication enabled: {ENABLE_AUTH}")

if ENABLE_AUTH:
    
    from authlib.integrations.flask_client import OAuth
    oauth = OAuth()  #Initialize OAuth only if auth is enabled
    

def create_flask_app(config_name=None):
    """ Create and configure the Flask app."""
    app = Flask(__name__, template_folder='../templates', static_folder='/app/dist/flask')
    # this was causing a number of issues with various routes, changing this here seems to be the best way to fix it
    # as there isn't a clear single point of front-end that would consistently guarantee fixing it
    app.url_map.strict_slashes = False
    app.after_request(add_safe_headers)
    
    app.config.update(
        SESSION_COOKIE_HTTPONLY=True,   # Prevent JavaScript from accessing cookies
        SESSION_COOKIE_SECURE=True,     # Only send cookies over HTTPS
        SESSION_COOKIE_SAMESITE="Lax"   # Prevent cross-site cookie usage
    )
    
    try:
        logger.info("Adding config.json details to app config")
        load_config(app, config_name, ENABLE_AUTH)
    except Exception as e:
        logger.exception(f"Error loading configuration: {e}")
        raise e
    
    try:
        logger.info("Creating base directory")
        create_base_directory(app)
    except Exception as e:
        logger.exception(f"Error creating base directory: {e}")
        raise e

    try:
        logger.info("Initializing app with db")
        db.init_app(app)
    except Exception as e:
        logger.exception(f"Error initializing database: {e}")
        raise e

    try:
        logger.info("Creating tables")
        with app.app_context():
            logger.info("Waiting for DB to set up")
            wait_for_database()
            if not tables_exist():
                logger.info("Creating database tables")
                db.create_all()
                logger.info("Created database tables")
            else:
                logger.info("Database tables already exist")

            if ENABLE_AUTH:
                try:
                    logger.info("Syncing users from Auth0 into the database...")
                    sync_auth0_users_to_db()

                    logger.info("Caching user-projects data...")
                    cache_user_projects()  # Cache the user-project mappings into Redis only when Auth is enabled

                except Exception as e:
                    logger.exception(f"Error during auth-related DB setup: {e}")
                    raise e
            
            # Routes registration and application setup
            logger.info("Registering the blueprint (register_app)")
            ProjectBlueprint.register_app(app)
    except OperationalError as oe:
        logger.exception(f"OperationalError: {oe}")
        raise oe
    except Exception as e:
        logger.exception(f"Error during app setup: {e}")
        raise e

    if ENABLE_AUTH:
        try:
            logger.info("Initializing OAuth for authentication")
            oauth.init_app(app)

            logger.info("Registering authentication before_request logic")
            register_before_request_auth(app)

            logger.info("Registering authentication routes")
            register_auth_routes(app)  # Register Auth0-related routes like /login and /callback

        except Exception as e:
            logger.exception(f"Error setting up authentication: {e}")
            raise e

    # Register other routes (base routes like /, /projects, etc.)
    try:
        # Register routes
        logger.info("Registering base routes: /, /projects, /create_project, /delete_project")
        register_routes(app, ENABLE_AUTH)
    except Exception as e:
        logger.exception(f"Error registering routes: {e}")
        raise e

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
            logger.info("Database is ready!")
            return
        except OperationalError as oe:
            logger.exception(f"OperationalError: {oe}. Database not ready, retrying in {delay} seconds... (Attempt {attempt + 1} of {max_retries})")
            time.sleep(delay)
        except Exception as e:
            logger.exception(f"An unexpected error occurred while waiting for the database: {e}")
            raise  # Re-raise the exception to be handled by the parent
            # ^^ should this be `raise e` instead?

    # If the loop completes without a successful connection
    error_message = "Error: Database did not become available in time."
    logger.error(error_message)
    raise TimeoutError(error_message)

# Load sensitive data from Docker secrets
def read_secret(secret_name):
    secret_path = f'/run/secrets/{secret_name}'
    try:
        with open(secret_path, 'r') as secret_file:
            return secret_file.read().strip()
    except FileNotFoundError as fnf_error:
        logger.exception(f"Error: Secret '{secret_name}' not found. {fnf_error}")
        raise  # Re-raise the exception to be handled by the parent.

def load_config(app, config_name=None, enable_auth=False):
    try:
        config_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.json')
        with open(config_file_path) as config_file:
            config = json.load(config_file)
            #app.config['PREFERRED_URL_SCHEME'] = 'https'
            app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = config.get('track_modifications', False)
            app.config['upload_folder'] = config.get('upload_folder', '')
            app.config['projects_base_dir'] = config.get('projects_base_dir', '')
            app.config['db_host'] = config.get('db_container', '')
            logger.info("Configuration loaded successfully!")
    except FileNotFoundError:
        logger.exception("Error: Configuration file not found.")
        raise
    except Exception as e:
        logger.exception(f"An unexpected error occurred: {e}")
        raise

    # Handle different environments
    try:
        if config_name == 'test':
            app.config['PREFERRED_URL_SCHEME'] = 'http'
            app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///:memory:'
        else:
            app.config['PREFERRED_URL_SCHEME'] = 'https'
            
            db_user = os.getenv('DB_USER') or read_secret('db_user')
            db_password = os.getenv('DB_PASSWORD') or read_secret('db_password')
            db_name = os.getenv('DB_NAME') or read_secret('db_name')
            db_host = os.getenv('DB_HOST') or app.config.get('db_host')

            if not all([db_user, db_password, db_name, db_host]):
                raise ValueError("Error: One or more required secrets or configurations are missing.")
            
            app.config['SQLALCHEMY_DATABASE_URI'] = f'postgresql://{db_user}:{db_password}@{db_host}/{db_name}'


            # Only configure Auth0 if ENABLE_AUTH is True
            if enable_auth:
                app.secret_key = os.getenv('FLASK_SECRET_KEY') or read_secret('flask_secret_key')
                auth0_domain = os.getenv('AUTH0_DOMAIN') or config.get('AUTH0_DOMAIN')
                auth0_client_id = os.getenv('AUTH0_CLIENT_ID') or config.get('AUTH0_CLIENT_ID')
                auth0_client_secret = os.getenv('AUTH0_CLIENT_SECRET') or read_secret("auth0_client_secret")

                if not all([auth0_domain, auth0_client_id, auth0_client_secret]):
                    raise ValueError("Error: Missing Auth0 configuration.")

                app.config['ENABLE_AUTH'] = True
                app.config['AUTH0_DOMAIN'] = auth0_domain
                app.config['AUTH0_CLIENT_ID'] = auth0_client_id
                app.config['AUTH0_CLIENT_SECRET'] = auth0_client_secret
                app.config['AUTH0_CALLBACK_URL'] = os.getenv('AUTH0_CALLBACK_URL') or config.get('AUTH0_CALLBACK_URL')
                app.config["AUTH0_PUBLIC_KEY_URI"] = os.getenv('AUTH0_PUBLIC_KEY_URI') or config.get('AUTH0_PUBLIC_KEY_URI')
                app.config["AUTH0_AUDIENCE"] = os.getenv('AUTH0_AUDIENCE') or config.get('AUTH0_AUDIENCE')
                app.config["AUTH0_DB_CONNECTION"] = os.getenv('AUTH0_DB_CONNECTION') or config.get('AUTH0_DB_CONNECTION')

                app.config["LOGIN_REDIRECT_URL"] = os.getenv('LOGIN_REDIRECT_URL') or config.get('LOGIN_REDIRECT_URL')

                app.config["SHIBBOLETH_LOGIN_URL"] = os.getenv('SHIBBOLETH_LOGIN_URL') or config.get('SHIBBOLETH_LOGIN_URL')
                app.config["SHIBBOLETH_LOGOUT_URL"] = os.getenv('SHIBBOLETH_LOGOUT_URL') or config.get('SHIBBOLETH_LOGOUT_URL')
            else:
                app.config['ENABLE_AUTH'] = False
    except Exception as e:
        logger.exception(f"An unexpected error occurred while configuring the database: {e}")
        raise  # Re-raise the exception to be handled by the parent
       
# Function to create base directory if it doesn't exist
def create_base_directory(app):
    try:
        base_dir = app.config.get('projects_base_dir', 'mdv')
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
            logger.info(f"Created base directory: {base_dir}")
        else:
            logger.info(f"Base directory already exists: {base_dir}")
    except Exception as e:
        logger.exception(f'Function create_base_directory Error: {e}')
        raise

def tables_exist():
        inspector = db.inspect(db.engine)
        #logger.info("printing table names")
        #print(inspector.get_table_names())
        return inspector.get_table_names()

def serve_projects_from_db(app):
    failed_projects: list[tuple[int, str | Exception]] = []
    try:
        # Get all projects from the database
        logger.info("Serving the projects present in both database and filesystem. Displaying the error if the path doesn't exist for a project")
        projects = Project.query.all()

        for project in projects:
            
            if project.is_deleted:
                logger.info(f"Project with ID {project.id} is soft deleted.")
                continue

            if os.path.exists(project.path):
                try:
                    p = MDVProject(dir=project.path, id=str(project.id), backend_db= True)
                    p.set_editable(True)
                    # todo: look up how **kwargs works and maybe have a shared app config we can pass around
                    p.serve(app=app, open_browser=False, backend_db=True)
                    logger.info(f"Serving project: {project.path}")

                    # Update or add files in the database to reflect the actual files in the filesystem
                    for root, dirs, files in os.walk(project.path):
                        for file_name in files:
                            full_file_path = os.path.join(root, file_name)

                            # Use the utility function to add or update the file in the database
                            try:
                                # Attempt to add or update the file in the database
                                FileService.add_or_update_file_in_project(
                                    file_name=file_name,
                                    file_path=full_file_path,
                                    project_id=project.id
                                )
                                #print(f"Processed file in DB: {file_name} at {full_file_path}")

                            except RuntimeError as file_error:
                                logger.exception(f"Failed to add or update file '{file_name}' in the database: {file_error}")


                except Exception as e:
                    logger.exception(f"Error serving project #{project.id}'{project.path}': {e}")
                    # don't `raise` here; continue serving other projects
                    # but keep track of failed projects & associated errors
                    # nb keeping track via project.id rather than instance of Project, because ORM seems to make that not work
                    failed_projects.append((project.id, e))
            else:
                e = f"Error serving project #{project.id}: path '{project.path}' does not exist."
                logger.error(e)
                failed_projects.append((project.id, e))
               
    except Exception as e:
        logger.exception(f"Error serving projects from database: {e}")
        raise
    logger.info(f"{len(failed_projects)} projects failed to serve. ({len(projects)} projects served successfully)")
    # nb using extend rather than replacing the list, but as of now I haven't made corresponding `serve_projects_from_filesytem` changes etc
    # so we really only expect this to run once, and the list to be empty
    ProjectService.failed_projects.extend(failed_projects)

def serve_projects_from_filesystem(app, base_dir):
    try:
        logger.info("Serving the projects present in filesystem but missing in database")
        logger.info(f"Scanning base directory: {base_dir}")

        # Get all project paths from the database
        projects_in_db = {project.path for project in Project.query.with_entities(Project.path).all()}
        logger.info(f"Project paths in DB: {projects_in_db}")

        # Get all project directories in the filesystem
        project_paths_in_fs = {os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))}
        logger.info(f"Project paths in filesystem: {project_paths_in_fs}")

        # Determine which project paths are in the filesystem but not in the database
        missing_project_paths = project_paths_in_fs - projects_in_db
        logger.info(f"Missing project paths: {missing_project_paths}")

        # Iterate over missing project paths to create and serve them
        for project_path in missing_project_paths:
            logger.info(f"Processing project path: {project_path}")
            
            if os.path.exists(project_path):
                try:
                    project_name = os.path.basename(project_path)

                    # Get the next ID from the database
                    next_id = db.session.query(db.func.max(Project.id)).scalar()
                    if next_id is None:
                        next_id = 1
                    else:
                        next_id += 1

                    p = MDVProject(dir=project_path, id= str(next_id), backend_db= True)
                    p.set_editable(True)
                    p.serve(app=app, open_browser=False, backend_db=True) 
                    logger.info(f"Serving project: {project_path}")

                    # Create a new Project record in the database with the default name
                    new_project = ProjectService.add_new_project(name=project_name, path=project_path)
                    if new_project is None:
                        raise ValueError(f"Failed to add project '{project_name}' to the database.")
                    
                    logger.info(f"Added project to DB: {new_project}")

                    # Rename directory to use project ID as folder name
                    project_id_str = str(new_project.id)
                    desired_path = os.path.join(app.config["projects_base_dir"], project_id_str)

                    if project_path != desired_path:
                        try:
                            # Rename the directory
                            os.rename(project_path, desired_path)
                            logger.info(f"Renamed project folder from {project_path} to {desired_path}")

                            # Update project path in DB
                            new_project.path = desired_path
                            db.session.commit()
                            logger.info(f"Updated project path in DB for project ID {new_project.id}")

                            # Also update local reference for downstream operations (like file sync)
                            project_path = desired_path

                        except Exception as rename_error:
                            logger.exception(f"Failed to rename project directory or update DB for project ID {new_project.id}: {rename_error}")

                    # Auth-related setup
                    if ENABLE_AUTH:
                        try:
                            sync_auth0_users_to_db()  # Sync users and assign permissions
                            logger.info("Synced Auth0 users after adding project.")

                            cache_user_projects() #update the cache
                        except Exception as auth_e:
                            logger.exception(f"Error syncing users or caching for {project_name}: {auth_e}")

                    
                    # Add files from the project directory to the database
                    for root, dirs, files in os.walk(project_path):
                        for file_name in files:
                            # Construct the full file path
                            full_file_path = os.path.join(root, file_name)
                            
                            # Use the full file path when adding or updating the file in the database
                            # Use the utility function to add or update the file in the database
                            try:
                                # Attempt to add or update the file in the database
                                FileService.add_or_update_file_in_project(
                                    file_name=file_name,
                                    file_path=full_file_path,
                                    project_id=new_project.id
                                )
                                #print(f"Processed file in DB: {file_name} at {full_file_path}")

                            except RuntimeError as file_error:
                                logger.exception(f"Failed to add or update file '{file_name}' in the database: {file_error}")
                except Exception as e:
                    logger.exception(f"In create_projects_from_filesystem: Error creating project at path '{project_path}': {e}")
                    raise
            else:
                logger.error(f"In create_projects_from_filesystem: Error - Project path '{project_path}' does not exist.")
                
    except Exception as e:
        logger.exception(f"In create_projects_from_filesystem: Error retrieving projects from database: {e}")
        raise


# Create the app object at the module level
app = create_flask_app()

with app.app_context():
    logger.info("Serving projects from database")
    serve_projects_from_db(app)
    logger.info("Starting - create_projects_from_filesystem")
    serve_projects_from_filesystem(app, app.config['projects_base_dir'])

if __name__ == '__main__':
    logger.info("Inside main..")
    #wait_for_database()
    logging.basicConfig(level=logging.INFO)

    app.run(host='0.0.0.0', debug=False, port=5055)