import os
import sys
import time
import json
import logging
from sqlalchemy import text, create_engine
from sqlalchemy.exc import OperationalError
from flask import Flask
from mdvtools.server import add_safe_headers
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprint
from mdvtools.dbutils.dbmodels import db, Project
from mdvtools.dbutils.routes import register_routes
from mdvtools.auth.register_auth_routes import register_auth_routes
from mdvtools.auth.authutils import register_before_request_auth, get_auth_provider, cache_user_projects
from mdvtools.dbutils.dbservice import ProjectService, FileService
from mdvtools.websocket import mdv_socketio
from mdvtools.logging_config import get_logger
from mdvtools.dbutils.server_options import get_server_options_for_db_projects
#this shouldn't be necessary in future
from psycogreen.gevent import patch_psycopg
patch_psycopg()

# Setup logging using the centralized logging module
logger = get_logger(__name__)

# Read environment flag for authentication
ENABLE_AUTH = os.getenv("ENABLE_AUTH", "0").lower() in ["1", "true", "yes"]
logger.info(f"Authentication enabled: {ENABLE_AUTH}")
oauth = None
if ENABLE_AUTH:
    
    from authlib.integrations.flask_client import OAuth
    oauth = OAuth()  #Initialize OAuth only if auth is enabled
    

def create_flask_app(config_name=None):
    """Create and configure the Flask app."""
    app = Flask(__name__, template_folder='../templates', static_folder='/app/dist/flask')

    if config_name == 'test':
        app.config['TESTING'] = True
        #Set a default for testing to avoid startup errors.
        app.config['DEFAULT_AUTH_METHOD'] = 'dummy'

    mdv_socketio(app)
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
            if config_name != 'test':
                wait_for_database(app)
            if not tables_exist():
                logger.info("Creating database tables")
                db.create_all()
                logger.info("Created database tables")
            else:
                logger.info("Database tables already exist")

            if ENABLE_AUTH:
                try:
                    logger.info("Syncing users from Auth provider into the database...")
                    auth_provider = get_auth_provider()
                    auth_provider.sync_users_to_db()

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
            assert(oauth is not None), "OAuth is not initialized"
            oauth.init_app(app)

            logger.info("Registering authentication before_request logic")
            register_before_request_auth(app)

            logger.info("Registering authentication routes")
            register_auth_routes(app)  # Register Auth0-related routes like /login and /callback

        except Exception as e:
            logger.exception(f"Error setting up authentication: {e}")
            raise e

    # Register other routes (base routes like /, /projects, /rescan_projects, etc.)
    # Note: Project management routes are now handled by ProjectManagerExtension
    try:
        # Register routes
        logger.info("Registering base routes: /, /projects, /rescan_projects")
        register_routes(app, ENABLE_AUTH)
    except Exception as e:
        logger.exception(f"Error registering routes: {e}")
        raise e

    # Register global routes from extensions
    try:
        logger.info("Registering global routes from extensions")
        from mdvtools.dbutils.server_options import get_server_options_for_db_projects
        options = get_server_options_for_db_projects(app)
        
        for extension in options.extensions:
            if hasattr(extension, 'register_global_routes'):
                logger.info(f"Registering global routes for extension: {extension.__class__.__name__}")
                extension.register_global_routes(app, app.config)
    except Exception as e:
        logger.exception(f"Error registering global routes from extensions: {e}")
        raise e

    return app


def wait_for_database(app):
    """Wait for the database to be ready before proceeding."""
    max_retries = 30
    delay = 5  # seconds
    db_name = os.getenv('DB_NAME')
    
    if not db_name:
        error_message = "Error: DB_NAME environment variable is not set"
        logger.error(error_message)
        raise ValueError(error_message)

    # First connect to postgres database to check/create our target database
    postgres_uri = 'postgresql://{}:{}@{}/postgres'.format(
        os.getenv('DB_USER'),
        os.getenv('DB_PASSWORD'),
        os.getenv('DB_HOST')
    )
    
    # Construct target URI once since environment variables don't change
    target_uri = 'postgresql://{}:{}@{}/{}'.format(
        os.getenv('DB_USER'),
        os.getenv('DB_PASSWORD'),
        os.getenv('DB_HOST'),
        db_name
    )
    
    # Update the database URI in the app config (this is the same as the initial URI from load_config)
    app.config['SQLALCHEMY_DATABASE_URI'] = target_uri
    
    for attempt in range(max_retries):
        try:
            # First connect to postgres database
            # db.create_engine does work, but not as well-typed, 
            # importing create_engine seems to be the standard way to do it
            # if not hasattr(db, 'create_engine'):
            #     raise Exception("db.create_engine not found")            
            # engine = db.create_engine(postgres_uri)
            engine = create_engine(postgres_uri, isolation_level="AUTOCOMMIT")
            connection = None
            try:
                connection = engine.connect()
                # Check if our target database exists
                result = connection.execute(
                    text("SELECT 1 FROM pg_database WHERE datname = :db_name"),
                    {"db_name": db_name}
                )
                exists = result.scalar()
                
                if not exists:
                    # Create the database if it doesn't exist
                    connection.execute(text(f'CREATE DATABASE "{db_name}"'))
                    logger.info(f"Created database: {db_name}")
                else:
                    logger.info(f"Database {db_name} already exists")
                if connection:
                    connection.close()
            finally:
                engine.dispose()
            
            # Test the connection to our target database
            # The engine should already be valid since the URI hasn't changed
            test_connection = None
            try:
                test_connection = db.engine.connect()
                test_connection.execute(text('SELECT 1'))
                logger.info("Successfully connected to target database")
            finally:
                if test_connection:
                    test_connection.close()
            
            logger.info("Database is ready!")
            return
            
        except OperationalError as oe:
            logger.exception(f"OperationalError: {oe}. Database not ready, retrying in {delay} seconds... (Attempt {attempt + 1} of {max_retries})")
            time.sleep(delay)
        except Exception as e:
            logger.exception(f"An unexpected error occurred while waiting for the database: {e}")
            raise

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
        # todo - as well as the default config stored in the repo, we can have a config.json for a given deployment
        config_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.json')
        with open(config_file_path) as config_file:
            config = json.load(config_file)
            #app.config['PREFERRED_URL_SCHEME'] = 'https'
            app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = config.get('track_modifications', False)
            app.config['upload_folder'] = config.get('upload_folder', '')
            app.config['projects_base_dir'] = config.get('projects_base_dir', '')
            app.config['db_host'] = config.get('db_container', '')
            # Control expensive per-file database sync during project serve
            sync_env = os.getenv('ENABLE_FILE_SYNC')
            if sync_env is not None:
                app.config['ENABLE_FILE_SYNC'] = sync_env.lower() in ["1", "true", "yes"]
            else:
                app.config['ENABLE_FILE_SYNC'] = config.get('enable_file_sync', False)
            # Allow extensions to be configured via user-provided JSON file for deployment flexibility
            # Check if external config path is provided via environment variable
            external_config_path = os.getenv('MDV_USER_CONFIG_PATH')
            if external_config_path:
                if os.path.exists(external_config_path):
                    try:
                        with open(external_config_path) as user_config_file:
                            user_config = json.load(user_config_file)
                            app.config['extensions'] = user_config.get('extensions', [])
                            logger.info(f"Loaded extensions from user config at {external_config_path}: {app.config['extensions']}")
                    except (json.JSONDecodeError, IOError) as e:
                        logger.warning(f"Could not load user config from {external_config_path}: {e}. Using default extensions.")
                        app.config['extensions'] = config.get('extensions', [])
                else:
                    logger.warning(f"User config file not found at {external_config_path}. Using default extensions.")
                    app.config['extensions'] = config.get('extensions', [])
            else:
                # No external config path specified, keep extensions empty
                app.config['extensions'] = []
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
                app.config['ENABLE_AUTH'] = True
                app.config["DEFAULT_AUTH_METHOD"] = os.getenv('DEFAULT_AUTH_METHOD') or config.get('DEFAULT_AUTH_METHOD')
                app.secret_key = os.getenv('FLASK_SECRET_KEY') or read_secret('flask_secret_key')
                
                # Check if the authentication method is 'auth0'
                if app.config["DEFAULT_AUTH_METHOD"] == "auth0":
                    auth0_domain = os.getenv('AUTH0_DOMAIN') or config.get('AUTH0_DOMAIN')
                    auth0_client_id = os.getenv('AUTH0_CLIENT_ID') or config.get('AUTH0_CLIENT_ID')
                    auth0_client_secret = os.getenv('AUTH0_CLIENT_SECRET') or read_secret("auth0_client_secret")

                    if not all([auth0_domain, auth0_client_id, auth0_client_secret]):
                        raise ValueError("Error: Missing Auth0 configuration.")

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

def is_valid_mdv_project(path: str):
    if not os.path.isdir(path):
        return False
    dir_list = os.listdir(path)
    return "views.json" in dir_list and "state.json" in dir_list and "datasources.json" in dir_list
        

def serve_projects_from_db(app):
    failed_projects: list[tuple[int, str | Exception]] = []
    try:
        # Get all projects from the database
        logger.info("Serving the projects present in both database and filesystem. Displaying the error if the path doesn't exist for a project")
        projects = Project.query.all()
        options = get_server_options_for_db_projects(app)

        for project in projects:
            
            if project.is_deleted:
                logger.info(f"Project with ID {project.id} is soft deleted.")
                continue

            if os.path.exists(project.path):
                try:
                    p = MDVProject(dir=project.path, id=str(project.id), backend_db= True)
                    # Sync DB access level from state.json if present, then serve with resulting editability
                    try:
                        state = p.state or {}
                        perm = (state.get('permission') or '').lower()
                        desired_level = 'editable' if perm == 'edit' else 'read-only' if perm == 'view' else None
                        if desired_level is not None and desired_level != getattr(project, 'access_level', None):
                            ProjectService.change_project_access(project.id, desired_level)
                            is_editable = (desired_level == 'editable')
                        else:
                            # Default to editable when access level is missing/unknown
                            is_editable = (project.access_level == 'editable') if getattr(project, 'access_level', None) else True
                        p.set_editable(is_editable)
                    except Exception:
                        # Favor editable by default on unexpected errors
                        p.set_editable(True)
                    # todo: look up how **kwargs works and maybe have a shared app config we can pass around
                    p.serve(options=options)
                    logger.info(f"Serving project: {project.path}")

                    # Optionally update/add files in DB to reflect filesystem
                    if app.config.get('ENABLE_FILE_SYNC', False):
                        for root, dirs, files in os.walk(project.path):
                            for file_name in files:
                                full_file_path = os.path.join(root, file_name)

                                try:
                                    FileService.add_or_update_file_in_project(
                                        file_name=file_name,
                                        file_path=full_file_path,
                                        project_id=project.id
                                    )
                                except RuntimeError as file_error:
                                    logger.exception(f"Failed to add or update file '{file_name}' in the database: {file_error}")
                    else:
                        logger.info("Skipping file sync for project %s (ENABLE_FILE_SYNC disabled)", project.id)


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
        project_paths_in_fs = {os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)) and d != 'lost+found'}
        logger.info(f"Project paths in filesystem: {project_paths_in_fs}")

        # Determine which project paths are in the filesystem but not in the database
        missing_project_paths = project_paths_in_fs - projects_in_db
        logger.info(f"Missing project paths: {missing_project_paths}")
        options = get_server_options_for_db_projects(app)

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
                    # Respect existing state.json permission if present; default to editable when unspecified
                    try:
                        state = p.state or {}
                        perm = (state.get('permission') or '').lower()
                        is_editable = True if perm == 'edit' else False if perm == 'view' else True
                        p.set_editable(is_editable)
                    except Exception:
                        p.set_editable(True)
                    p.serve(options=options) 
                    logger.info(f"Serving project: {project_path}")

                    # Create a new Project record in the database with the default name
                    new_project = ProjectService.add_new_project(name=project_name, path=project_path)
                    if new_project is None:
                        raise ValueError(f"Failed to add project '{project_name}' to the database.")
                    
                    logger.info(f"Added project to DB: {new_project}")
                    # One-time sync: initialize DB access_level from state.json.permission
                    try:
                        state = p.state or {}
                        perm = (state.get('permission') or '').lower()
                        desired_level = 'editable' if perm == 'edit' else 'read-only' if perm == 'view' else None
                        if desired_level is not None:
                            ProjectService.change_project_access(new_project.id, desired_level)
                    except Exception:
                        pass

                    # Rename directory to use project ID as folder name
                    """
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
                    """
                    # Auth-related setup
                    if ENABLE_AUTH:
                        try:
                            auth_provider = get_auth_provider()
                            auth_provider.sync_users_to_db()  # Sync users and assign permissions
                            logger.info("Synced Auth users after adding project.")

                            cache_user_projects() #update the cache
                        except Exception as auth_e:
                            logger.exception(f"Error syncing users or caching for {project_name}: {auth_e}")

                    
                    # Optionally add files from the project directory to the database
                    if app.config.get('ENABLE_FILE_SYNC', False):
                        for root, dirs, files in os.walk(project_path):
                            for file_name in files:
                                full_file_path = os.path.join(root, file_name)
                                try:
                                    FileService.add_or_update_file_in_project(
                                        file_name=file_name,
                                        file_path=full_file_path,
                                        project_id=new_project.id
                                    )
                                except RuntimeError as file_error:
                                    logger.exception(f"Failed to add or update file '{file_name}' in the database: {file_error}")
                    else:
                        logger.info("Skipping file sync for new project %s (ENABLE_FILE_SYNC disabled)", new_project.id)
                except Exception as e:
                    logger.exception(f"In create_projects_from_filesystem: Error creating project at path '{project_path}': {e}")
                    raise
            else:
                logger.error(f"In create_projects_from_filesystem: Error - Project path '{project_path}' does not exist.")
                
    except Exception as e:
        logger.exception(f"In create_projects_from_filesystem: Error retrieving projects from database: {e}")
        raise


# Create the app object at the module level
try:
    app = create_flask_app()
    
    with app.app_context():
        logger.info("Serving projects from database")
        serve_projects_from_db(app)
        logger.info("Starting - create_projects_from_filesystem")
        serve_projects_from_filesystem(app, app.config['projects_base_dir'])
except Exception as e:
    logger.exception(f"Error during app initialization: {e}")
    app = None

if __name__ == '__main__':
    logger.info("Inside main..")
    #wait_for_database()

    if app is None:
        logger.error("App initialization failed, cannot start server")
        sys.exit(1)

    app.run(host='0.0.0.0', debug=False, port=5055)