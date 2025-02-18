import os
import time
import json
import shutil
from sqlalchemy import text
from sqlalchemy.exc import OperationalError
from flask import Flask, render_template, jsonify, request
from flask import Flask, render_template, jsonify, request
#from flask_sqlalchemy import SQLAlchemy
# import threading
# from flask import Flask, render_template, jsonify, request
from mdvtools.server import add_safe_headers
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprint_v2 as ProjectBlueprint
from mdvtools.dbutils.dbmodels import db, Project
#from mdvtools.dbutils.routes import register_global_routes
from mdvtools.dbutils.dbservice import ProjectService, FileService
from flask import redirect, url_for, session, jsonify

# Read environment flag for authentication
ENABLE_AUTH = os.getenv("ENABLE_AUTH", "0").lower() in ["1", "true", "yes"]

if ENABLE_AUTH:
    
        from authlib.integrations.flask_client import OAuth
        from mdvtools.auth.auth0_provider import Auth0Provider

        oauth = OAuth()  # Initialize OAuth only if auth is enabled
    #except ImportError:
    #    print("Auth library not found. Ensure poetry installs `auth` dependencies when ENABLE_AUTH=1.")
    #    exit(1)  # Fail early if auth is enabled but libraries are missing


def create_flask_app(config_name=None):
    """Create and configure the Flask app."""
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
    app.secret_key = os.getenv('FLASK_SECRET_KEY', 'a-very-secret-key')

    
    try:
        print("** Adding config.json details to app config")
        load_config(app, config_name, ENABLE_AUTH)
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
                print("************** Created database tables")
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

    # Register OAuth with the app
    if ENABLE_AUTH:
        try:
            print("Initializing OAuth for authentication")
            oauth.init_app(app)
        except Exception as e:
            print(f"Error initializing OAuth: {e}")
            exit(1)

    # Global variable for Auth0 provider
    auth0_provider = None  # This should be set when Auth0 routes are registered

    def is_authenticated():
        """Check if the user is authenticated, considering Auth0 and Shibboleth."""
        if not ENABLE_AUTH:
            return True  # Authentication is disabled, allow access
        
        if session.get("auth_method") == "shibboleth":
            return True  # Shibboleth users are already authenticated
        
        return "token" in session and session["token"]

    # Authentication check function
    def is_authenticated_token():
        """Check if the user is authenticated (works for both Auth0 and Shibboleth)."""
        # Check if auth0_provider is initialized
        if not auth0_provider:
            print("<<<<<<1")
            # If auth0_provider is not available, return False as we can't check authentication for Auth0 yet
            return False
        
        # Check if Auth0 is enabled and if the token is in session
        if session.get('auth_method') == 'auth0' and 'token' in session:
            print("<<<<<<2")
            # For Auth0, we need to validate the token to ensure it's not expired or invalid
            try:
                # Validate the token with Auth0 (this will depend on your token structure and verification method)
                if auth0_provider.is_token_valid(session['token']):
                    print("--------88888")
                    return True  # Token is valid, user is authenticated
                else:
                    # Token is invalid or expired
                    print("Auth0 token is invalid or expired.")
                    session.clear()  # Clear session if token is invalid
                    return False
            except Exception as e:
                print(f"Error during Auth0 token validation: {e}")
                session.clear()  # Clear session if there was an error during validation
                return False
        
        # Check for Shibboleth authentication
        elif session.get('auth_method') == 'shibboleth':
            # For Shibboleth, if the token is not available in the session, we cannot validate like Auth0
            # So, we assume that the presence of 'auth_method' is the sign of successful authentication
            # You may adjust this depending on your specific Shibboleth setup.
            if 'auth_method' in session and session['auth_method'] == 'shibboleth':
                return True  # Assume user is authenticated after successful Shibboleth login
            
            return False  # If no 'auth_method' is found for Shibboleth, user is not authenticated
        
        return False  # Default to False if neither Auth0 nor Shibboleth is used


    # Whitelist of routes that do not require authentication
    whitelist_routes = [
        '/login_dev',
        '/login_sso',
        '/login',
        '/callback',
        '/favicon.ico',        # Allow access to favicon
        '/flask/js/',  # Allow access to login JS
        '/static',
        '/flask/assets',
        '/flask/img'
    ]

    @app.before_request
    def enforce_authentication():
        """Redirect unauthenticated users to login if required."""
        if not ENABLE_AUTH:
            print(":::::1")
            return None  # Skip authentication check if auth is disabled

        requested_path = request.path
        if any(requested_path.startswith(route) for route in whitelist_routes):
            print(":::::2")
            return None  # Allow access to whitelisted routes

        if not is_authenticated():
            print(":::::3")
            redirect_uri = app.config["LOGIN_REDIRECT_URL"]
            print(f"Unauthorized access attempt to {requested_path}. Redirecting to /login_dev.")
            return redirect(redirect_uri)

        return None

    if ENABLE_AUTH:
        try:
            print("Registering authentication routes")
            register_auth0_routes(app)  # Register Auth0-related routes like /login and /callback
        except Exception as e:
            print(f"Error registering authentication routes: {e}")
            exit(1)

    # Register other routes (base routes like /, /projects, etc.)
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
            print(f"OperationalError: {oe}. Database not ready, retrying in {delay} seconds... (Attempt {attempt + 1} of {max_retries})")
            time.sleep(delay)
        except Exception as e:
            print(f"An unexpected error occurred while waiting for the database: {e}")
            raise  # Re-raise the exception to be handled by the parent
            # ^^ should this be `raise e` instead?

    # If the loop completes without a successful connection
    error_message = "Error: Database did not become available in time."
    print(error_message)
    raise TimeoutError(error_message)



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
            app.config['PREFERRED_URL_SCHEME'] = 'https'
            
            # Load sensitive data from Docker secrets
            def read_secret(secret_name):
                secret_path = f'/run/secrets/{secret_name}'
                try:
                    with open(secret_path, 'r') as secret_file:
                        return secret_file.read().strip()
                except FileNotFoundError as fnf_error:
                    print(f"Error: Secret '{secret_name}' not found. {fnf_error}")
                    raise  # Re-raise the exception to be handled by the parent

            db_user = os.getenv('DB_USER') or read_secret('db_user')
            db_password = os.getenv('DB_PASSWORD') or read_secret('db_password')
            db_name = os.getenv('DB_NAME') or read_secret('db_name')
            db_host = os.getenv('DB_HOST') or app.config.get('db_host')

            if not all([db_user, db_password, db_name, db_host]):
                raise ValueError("Error: One or more required secrets or configurations are missing.")
            
            app.config['SQLALCHEMY_DATABASE_URI'] = f'postgresql://{db_user}:{db_password}@{db_host}/{db_name}'


            # Only configure Auth0 if ENABLE_AUTH is True
            if enable_auth:
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

                app.config["LOGIN_REDIRECT_URL"] = os.getenv('LOGIN_REDIRECT_URL') or config.get('LOGIN_REDIRECT_URL')

                app.config["SHIBBOLETH_LOGIN_URL"] = os.getenv('SHIBBOLETH_LOGIN_URL') or config.get('SHIBBOLETH_LOGIN_URL')
                app.config["SHIBBOLETH_LOGOUT_URL"] = os.getenv('SHIBBOLETH_LOGOUT_URL') or config.get('SHIBBOLETH_LOGOUT_URL')

                


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
    failed_projects: list[tuple[int, str | Exception]] = []
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
                    p = MDVProject(dir=project.path, id=str(project.id), backend_db= True)
                    p.set_editable(True)
                    # todo: look up how **kwargs works and maybe have a shared app config we can pass around
                    p.serve(app=app, open_browser=False, backend_db=True)
                    print(f"Serving project: {project.path}")

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
                                print(f"Failed to add or update file '{file_name}' in the database: {file_error}")


                except Exception as e:
                    print(f"Error serving project #{project.id}'{project.path}': {e}")
                    # don't `raise` here; continue serving other projects
                    # but keep track of failed projects & associated errors
                    # nb keeping track via project.id rather than instance of Project, because ORM seems to make that not work
                    failed_projects.append((project.id, e))
            else:
                e = f"Error serving project #{project.id}: path '{project.path}' does not exist."
                print(e)
                failed_projects.append((project.id, e))
               
    except Exception as e:
        print(f"Error serving projects from database: {e}")
        raise
    print(f"{len(failed_projects)} projects failed to serve. ({len(projects)} projects served successfully)")
    # nb using extend rather than replacing the list, but as of now I haven't made corresponding `serve_projects_from_filesytem` changes etc
    # so we really only expect this to run once, and the list to be empty
    ProjectService.failed_projects.extend(failed_projects)

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

                    p = MDVProject(dir=project_path, id= str(next_id), backend_db= True)
                    p.set_editable(True)
                    p.serve(app=app, open_browser=False, backend_db=True) 
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
                                print(f"Failed to add or update file '{file_name}' in the database: {file_error}")
                except Exception as e:
                    print(f"In create_projects_from_filesystem: Error creating project at path '{project_path}': {e}")
                    raise
            else:
                print(f"In create_projects_from_filesystem: Error - Project path '{project_path}' does not exist.")
                
    except Exception as e:
        print(f"In create_projects_from_filesystem: Error retrieving projects from database: {e}")
        raise


# The function that registers the Auth0 routes
def register_auth0_routes(app):
    """
    Registers the Auth0 routes like login, callback, logout, etc. to the Flask app,
    with centralized and route-specific error handling.
    """
    print("Registering AUTH routes...")

    try:
        # Initialize the Auth0Provider
        auth0_provider = Auth0Provider(
            app,
            oauth=oauth,
            client_id=app.config['AUTH0_CLIENT_ID'],
            client_secret=app.config['AUTH0_CLIENT_SECRET'],
            domain=app.config['AUTH0_DOMAIN']
        )

        # Route for login (redirects to Auth0 for authentication)
        @app.route('/login')
        def login():
            try:
                print("$$$$$$$$$$$$$$$ app-login")
                session.clear()  
                return auth0_provider.login()
            except Exception as e:
                print(f"In register_auth0_routes : Error during login: {e}")
                return jsonify({"error": "Failed to start login process."}), 500

        # Route for the callback after login (handles the callback from Auth0)
        @app.route('/callback')
        def callback():
            try:
                print("$$$$$$$$$$$$$$$ app-callback")
                code = request.args.get('code')  # Get the code from the callback URL
                if not code:
                    print("Missing 'code' parameter in the callback URL.")
                    session.clear()  # Clear session if there's no code
                    return jsonify({"error": "Authorization code not provided."}), 400
                
                print("$$$$$$$$$$$$$$$ app-callback  1")
                access_token = auth0_provider.handle_callback()
                if not access_token:  # If token retrieval fails, prevent redirecting
                    print("Authentication failed: No valid token received.")
                    session.clear()  # Clear session in case of failure
                    return jsonify({"error": "Authentication failed."}), 401
                
                print(" $$$$$$$$$$$$$$$ app-callback 2")
                return redirect(url_for('index'))  # Redirect to the home page or any protected page
            except Exception as e:
                print(f"In register_auth0_routes : Error during callback: {e}")
                session.clear()  # Clear session on error
                return jsonify({"error": "Failed to complete authentication process."}), 500

        # Route for logout (clears the session and redirects to home)
        @app.route('/logout')
        def logout():
            try:
                # Check what authentication method was used (Auth0 or Shibboleth)
                auth_method = session.get('auth_method', None)

                if auth_method == 'auth0':
                    # If the user logged in via Auth0, log them out from Auth0
                    return auth0_provider.logout()
                    

                # If the user logged in via Shibboleth, redirect to Shibboleth IdP's logout URL
                elif auth_method == 'shibboleth':
                    # Shibboleth does not handle the session clearing, so we first clear the session
                    session.clear()
                    # Then, redirect to the Shibboleth IdP logout URL
                    shibboleth_logout_url = app.config.get('SHIBBOLETH_LOGOUT_URL', None)

                    if shibboleth_logout_url:
                        # Redirect to the provided Shibboleth IdP logout URL
                        return redirect(shibboleth_logout_url)
                    else:
                        # If no Shibboleth logout URL is configured, return an error
                        return jsonify({"error": "Shibboleth logout URL not provided."}), 500

                # Clear the session data after logging out from either Auth0 or Shibboleth
                session.clear()

                # No need to redirect here if auth0_provider.logout() already handles redirection
                return jsonify({"message": "Logged out successfully"}), 200

            except Exception as e:
                print(f"In register_auth0_routes: Error during logout: {e}")
                session.clear()
                return jsonify({"error": "Failed to log out."}), 500


        # You can also add a sample route to check the user's profile or token
        @app.route('/profile')
        def profile():
            try:
                token = session.get('token')
                if token:
                    user_info = auth0_provider.get_user(token)
                    return jsonify(user_info)
                else:
                    print("Token not found in session.")
                    return jsonify({"error": "Not authenticated."}), 401
            except Exception as e:
                print(f"In register_auth0_routes: Error during profile retrieval: {e}")
                return jsonify({"error": "Failed to retrieve user profile."}), 500
        

        @app.route('/login_sso')
        def login_sso():
            """Redirect user to Shibboleth-protected login page on Apache."""
            try:
                # Clear any existing session data to ensure we start with a fresh session
                session.clear()
                
                # Store the authentication method as Shibboleth
                session["auth_method"] = "shibboleth"  # Indicate Shibboleth login

                # Check if the Shibboleth login URL is provided in the environment
                shibboleth_login_url = app.config.get('SHIBBOLETH_LOGIN_URL', None)

                if shibboleth_login_url:
                    # Redirect the user to Shibboleth login page if the URL is configured
                    print("Redirecting to Shibboleth login page...")
                    return redirect(shibboleth_login_url)
                else:
                    # If Shibboleth URL is not provided, inform the user
                    print("Shibboleth login URL not provided.")
                    return jsonify({"error": "Shibboleth login URL not provided."}), 500

            except Exception as e:
                # In case of error, clear the session and handle the error
                session.clear()  # Ensure session is cleared in case of failure
                print(f"In login_sso: Error during login: {e}")
                return jsonify({"error": "Failed to start login process using SSO."}), 500

        print("Auth0 routes registered successfully!")

    except Exception as e:
        print(f"Error registering AUTH routes: {e}")
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

        @app.route('/login_dev')
        def login_dev():
            return render_template('login.html')
        print("Route registered: /login_dev")

        @app.route('/projects')
        def get_projects():
            print('/projects queried...')
            try:
                # Query the database to get all projects that aren't deleted
                projects = ProjectService.get_active_projects()
                
                # Format each project with its id, name, and last modified timestamp as a string
                project_list = [
                    {
                        "id": p.id,
                        "name": p.name,
                        "lastModified": p.update_timestamp.strftime('%Y-%m-%d %H:%M:%S')  # Format datetime as string
                    }
                    for p in projects
                ]
                # Return the list of projects as JSON
                return jsonify(project_list)
            
            except Exception as e:
                print(f"In register_routes - /projects : Error retrieving projects: {e}")
                return jsonify({"status": "error", "message": str(e)}), 500

        print("Route registered: /projects")

        @app.route("/create_project", methods=["POST"])
        def create_project():
            project_path = None
            next_id = None
            try:
                print("Creating project")
                
                # Get the next available ID
                next_id = ProjectService.get_next_project_id()
                if next_id is None:
                    print("In register_routes: Error- Failed to determine next project ID from db")
                    return jsonify({"status": "error", "message": "Failed to determine next project ID from db"}), 500

                # Create the project directory path
                project_path = os.path.join(app.config['projects_base_dir'], str(next_id))

                # Create and serve the MDVProject
                try:
                    print("Creating and serving the new project")
                    p = MDVProject(project_path, backend_db= True)
                    p.set_editable(True)
                    p.serve(app=app, open_browser=False, backend_db=True)
                except Exception as e:
                    print(f"In register_routes: Error serving MDVProject: {e}")
                    return jsonify({"status": "error", "message": "Failed to serve MDVProject"}), 500

                # Create a new Project record in the database with the path
                print("Adding new project to the database")
                new_project = ProjectService.add_new_project(path=project_path)

                if new_project:
                    return jsonify({"id": new_project.id, "name": new_project.name, "status": "success"})
                
                

            except Exception as e:
                print(f"In register_routes - /create_project : Error creating project: {e}")
                print("started rollabck")
                # Rollback: Clean up the projects filesystem directory if it was created
                if project_path and os.path.exists(project_path):
                    try:
                        shutil.rmtree(project_path)
                        print("In register_routes -/create_project : Rolled back project directory creation as db entry is not added")
                    except Exception as cleanup_error:
                        print(f"In register_routes -/create_project : Error during rollback cleanup: {cleanup_error}")

                # Optional: Remove project routes from Flask app if needed
                if next_id is not None and str(next_id) in ProjectBlueprint.blueprints:
                    del ProjectBlueprint.blueprints[str(next_id)]
                    print("In register_routes -/create_project : Rolled back ProjectBlueprint.blueprints as db entry is not added")
                
                return jsonify({"status": "error", "message": str(e)}), 500

        print("Route registered: /create_project")

        @app.route("/delete_project/<project_id>", methods=["DELETE"])
        def delete_project(project_id: int):
            #project_removed_from_blueprints = False
            try:
                print(f"Deleting project '{project_id}'")
                
                # Find the project by ID 
                project = ProjectService.get_project_by_id(project_id)

                if project is None:
                    print(f"In register_routes - /delete_project Error: Project with ID {project_id} not found in database")
                    return jsonify({"status": "error", "message": f"Project with ID {project_id} not found in database"}), 404

                
                # Check if the project is editable before attempting to delete
                if project.access_level != 'editable':
                    print(f"In register_routes - /delete_project Error: Project with ID {project_id} is not editable.")
                    return jsonify({"status": "error", "message": "This project is not editable and cannot be deleted."}), 403

                # Remove the project from the ProjectBlueprint.blueprints dictionary
                if str(project_id) in ProjectBlueprint.blueprints:
                    del ProjectBlueprint.blueprints[str(project_id)]
                    #project_removed_from_blueprints = True  # Mark as removed
                    print(f"In register_routes - /delete_project : Removed project '{project_id}' from ProjectBlueprint.blueprints")
                
                # Soft delete the project
                delete_status = ProjectService.soft_delete_project(project_id)

                if delete_status:
                    return jsonify({"status": "success"})
                else:
                    print("In register_routes - /delete_project Error: Failed to soft delete project in db")
                    return jsonify({"status": "error", "message": "Failed to soft delete project in db"}), 500

            except Exception as e:
                print(f"In register_routes - /delete_project: Error deleting project '{project_id}': {e}")
                return jsonify({"status": "error", "message": str(e)}), 500

        print("Route registered: /delete_project/<project_id>")

        @app.route("/projects/<int:project_id>/rename", methods=["PUT"])
        def rename_project(project_id: int):
            # Retrieve the new project name from the multipart/form-data payload
            new_name = request.form.get("name")
            
            if not new_name:
                return jsonify({"status": "error", "message": "New name not provided"}), 400
            
            try:
                # Check if the project exists
                project = ProjectService.get_project_by_id(project_id)
                if project is None:
                    print(f"In register_routes - /rename_project Error: Project with ID {project_id} not found in database")
                    return jsonify({"status": "error", "message": f"Project with ID {project_id} not found in database"}), 404
                
                # Check if the project is editable before attempting to rename
                if project.access_level != 'editable':
                    print(f"In register_routes - /rename_project Error: Project with ID {project_id} is not editable.")
                    return jsonify({"status": "error", "message": "This project is not editable and cannot be renamed."}), 403
      
                # Attempt to rename the project
                rename_status = ProjectService.update_project_name(project_id, new_name)

                if rename_status:
                    return jsonify({"status": "success", "id": project_id, "new_name": new_name}), 200
                else:
                    print(f"In register_routes - /rename_project Error: The project with ID '{project_id}' not found in db")
                    return jsonify({"status": "error", "message": f"Failed to rename project '{project_id}' in db"}), 500

            except Exception as e:
                print(f"In register_routes - /rename_project : Error renaming project '{project_id}': {e}")
                return jsonify({"status": "error", "message": str(e)}), 500

        print("Route registered: /projects/<int:project_id>/rename")

        @app.route("/projects/<int:project_id>/access", methods=["PUT"])
        def change_project_access(project_id):
            """API endpoint to change the access level of a project."""
            try:
                # Get the new access level from the request
                new_access_level = request.form.get("type")

                # Validate the new access level
                if new_access_level not in ["read-only", "editable"]:
                    return jsonify({"status": "error", "message": "Invalid access level. Must be 'read-only' or 'editable'."}), 400

                # Call the service method to change the access level
                access_level, message, status_code = ProjectService.change_project_access(project_id, new_access_level)

                if access_level is None:
                    return jsonify({"status": "error", "message": message}), status_code

                return jsonify({"status": "success", "access_level": access_level}), 200

            except Exception as e:
                print(f"In register_routes - /access : Unexpected error while changing access level for project '{project_id}': {e}")
                return jsonify({"status": "error", "message": "An unexpected error occurred."}), 500
        
        print("Route registered: /projects/<int:project_id>/access")

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