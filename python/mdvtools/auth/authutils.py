import logging
from functools import wraps
from flask import jsonify, session, request,redirect, current_app

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# in_memory_cache.py
user_cache = {}  # key: auth0_id -> user details
user_project_cache = {}  # key: user_id -> project permissions
all_users_cache = []  # list of all user summaries
active_projects_cache = []

def is_authenticated():
    """Validate the current user via Auth0 or Shibboleth."""
    ENABLE_AUTH = current_app.config.get("ENABLE_AUTH", False)
    
    if not ENABLE_AUTH:
        return True

    auth_method = session.get("auth_method")
    if auth_method == "auth0":
        user, error_response = validate_and_get_user()
    elif auth_method == "shibboleth":
        user, error_response = validate_sso_user(request)
    else:
        return False

    return user is not None and error_response is None

def register_before_request_auth(app):
    """Attach the before_request auth logic to the Flask app."""

    whitelist_routes = [
        '/login_dev',
        '/login_sso',
        '/login',
        '/callback',
        '/favicon.ico',
        '/flask/js/',
        '/static',
        '/flask/assets',
        '/flask/img'
    ]

    @app.before_request
    def enforce_authentication():
        ENABLE_AUTH = app.config.get("ENABLE_AUTH", False)

        if not ENABLE_AUTH:
            return None

        requested_path = request.path
        if any(requested_path.startswith(route) for route in whitelist_routes):
            return None

        if not is_authenticated():
            redirect_uri = app.config.get("LOGIN_REDIRECT_URL", "/login_dev")
            logger.info(f"Unauthorized access to {requested_path}. Redirecting.")
            return redirect(redirect_uri)

        return None

def sync_auth0_users_to_db():
    from mdvtools.dbutils.dbservice import UserService, UserProjectService
    from mdvtools.dbutils.dbmodels import db, Project
    from auth0.management import Auth0
    from auth0.authentication import GetToken

    """
    Syncs users from Auth0 to the application's database using UserService and UserProjectService.
    """
    try:
        # Load Auth0 config from app
        auth0_domain = current_app.config['AUTH0_DOMAIN']
        client_id = current_app.config['AUTH0_CLIENT_ID']
        client_secret = current_app.config['AUTH0_CLIENT_SECRET']
        auth0_db_connection = current_app.config['AUTH0_DB_CONNECTION']
        audience = f"https://{auth0_domain}/api/v2/"

        # Get Auth0 Management API token
        get_token = GetToken(domain=auth0_domain, client_id=client_id, client_secret=client_secret)
        mgmt_api_token = get_token.client_credentials(audience=audience)["access_token"]
        auth0 = Auth0(auth0_domain, mgmt_api_token)

        # Fetch users from Auth0 connection
        users = auth0.users.list(q=f'identities.connection:"{auth0_db_connection}"', per_page=100)

        for user in users['users']:
            email = user.get('email', '')
            auth0_id = user['user_id']

            # Use UserService to add or update user
            db_user = UserService.add_or_update_user(
                email=email,
                auth0_id=auth0_id
            )

            # Fetch user's roles to determine admin status
            roles = auth0.users.list_roles(auth0_id)
            is_admin = any(role['name'] == 'admin' for role in roles['roles'])

            # Update admin status
            db_user.is_admin = is_admin
            db.session.commit()

            if is_admin:
                # Assign all projects to this user as owner via UserProjectService
                for project in Project.query.all():
                    UserProjectService.add_or_update_user_project(
                        user_id=db_user.id,
                        project_id=project.id,
                        is_owner=True
                    )

        logger.info(f"Synced users from Auth0.")

    except Exception as e:
        logger.exception(f"sync_auth0_users_to_db: An unexpected error occurred: {e}")
        raise

def cache_user_projects():
    from mdvtools.dbutils.dbmodels import User, UserProject
    from mdvtools.dbutils.dbservice import ProjectService
    """
    Caches user details and their associated project permissions in memory.
    """
    try:
        logger.info("Caching user details and project permissions...")

        users = User.query.all()
        all_users_cache.clear()
        user_cache.clear()
        user_project_cache.clear()
        active_projects_cache.clear()

        for user in users:
            user_data = {
                "id": user.id,
                "auth0_id": user.auth0_id,
                "email": user.email,
                "is_admin": user.is_admin
            }
            user_cache[user.auth0_id] = user_data
            all_users_cache.append(user_data)

            user_projects = UserProject.query.filter_by(user_id=user.id).all()
            project_permissions = {
                up.project_id: {
                    "can_read": up.can_read,
                    "can_write": up.can_write,
                    "is_owner": up.is_owner
                }
                for up in user_projects
            }
            user_project_cache[user.id] = project_permissions

        active_projects = ProjectService.get_active_projects()
        active_projects_cache[:] = active_projects

        logger.info(f"Cached users and their project permissions in memory.")
        return True

    except Exception as e:
        logger.exception(f"Error caching user projects: {e}")
        return False
    
def update_cache(user_id=None, project_id=None, user_data=None, project_data=None, permissions=None):
    """
    Updates the in-memory caches (user_cache, user_project_cache, all_users_cache, active_projects_cache)
    for the provided user and project details, including any changes in permissions.

    :param user_id: The ID of the user whose cache needs to be updated.
    :param project_id: The ID of the project to be added or updated in the cache.
    :param user_data: A dictionary containing the user's details (if updating the user cache).
    :param project_data: A dictionary containing the project's details (if updating the project cache).
    :param permissions: A dictionary containing the user's permissions for a project (optional).
    """
    try:
        # Step 1: Update user_cache (user details)
        if user_data and user_id:
            user_cache[user_id] = user_data  # Update or add the user in the cache
            
            # Ensure the user is in all_users_cache if they are not already
            if user_data not in all_users_cache:
                all_users_cache.append(user_data)
        
        # Step 2: Update user_project_cache (user-project relationship)
        if user_id and project_id and permissions is not None:
            # Create or update user project permissions in cache
            if user_id not in user_project_cache:
                user_project_cache[user_id] = {}

            user_project_cache[user_id][project_id] = permissions  # Update or add project permissions
            
        # Step 3: Update active_projects_cache (active project details)
        if project_data and project_id:
            # Search for the project in the cache
            existing_project = next((p for p in active_projects_cache if p["id"] == project_id), None)

            if existing_project:
                # Project exists in cache, update only the changed fields
                existing_project["name"] = project_data.get("name", existing_project["name"])
                existing_project["lastModified"] = project_data.get("lastModified", existing_project["lastModified"])
                existing_project["thumbnail"] = project_data.get("thumbnail", existing_project["thumbnail"])
                logger.info(f"Updated project {project_id} in active projects cache.")
            else:
                # Project does not exist in cache, append a new entry
                project_entry = {
                    "id": project_data["id"],
                    "name": project_data["name"],
                    "lastModified": project_data["lastModified"],
                    "thumbnail": project_data["thumbnail"]
                }
                active_projects_cache.append(project_entry)
                logger.info(f"Added new project {project_id} to active projects cache.")

        logger.info("Cache successfully updated.")
    
    except Exception as e:
        logger.exception(f"Error updating cache: {e}")

def validate_and_get_user():
    from mdvtools.auth.auth0_provider import Auth0Provider
    from mdvtools.dbutils.mdv_server_app import oauth
    from mdvtools.dbutils.dbmodels import User
    """
    Validates the Auth0 token from the session and retrieves the user from cache.
    
    :param app: Flask app instance
    :return: Tuple (user object as dict, error response)
    """
    try:
        # Check if user information is already cached in session
        if 'user' in session:
            return session['user'], None  # Return the user from the session cache

        # Initialize Auth0 provider
        auth0_provider = Auth0Provider(
            current_app,
            oauth=oauth,
            client_id=current_app.config['AUTH0_CLIENT_ID'],
            client_secret=current_app.config['AUTH0_CLIENT_SECRET'],
            domain=current_app.config['AUTH0_DOMAIN']
        )

        # Retrieve the token from session
        token = auth0_provider.get_token()

        if not token:
            return None, (jsonify({"error": "Authentication required"}), 401)

        # Validate token
        if not auth0_provider.is_token_valid(token):
            return None, (jsonify({"error": "Invalid or expired token"}), 401)

        # Retrieve user info from Auth0
        user_info = auth0_provider.get_user({"access_token": token})

        if not user_info:
            return None, (jsonify({"error": "User not found"}), 404)

        # Get Auth0 user ID
        auth0_id = user_info.get("sub")

        # Fetch user directly from in-memory cache
        cached_user = user_cache.get(auth0_id)
        if cached_user:
            session['user'] = cached_user
            return cached_user, None  # Return cached user

        # This should rarely happen since all users were cached at app startup.
        # But if the user is missing from Redis, fallback to DB and log the issue.
        logger.info(f"User {auth0_id} not found in cache, falling back to DB!")

        # Query the user from the database if not in cache
        user = User.query.filter_by(auth0_id=auth0_id).first()
        if not user:
            return None, (jsonify({"error": "User not found"}), 404)

        # Add the user to the in-memory cache
        user_data = {"id": user.id, "auth0_id": user.auth0_id, "email": user.email, "is_admin": user.is_admin}
        user_cache[auth0_id] = user_data

        # Cache the user data in session for future use
        session['user'] = user_data

        return user_data, None

    except Exception as e:
        logger.exception(f"Error during authentication: {e}")
        return None, (jsonify({"error": "Internal server error - user not validated"}), 500)
    

def validate_sso_user(request):
    """
    Validates the SSO-authenticated user by checking required headers.
    Returns the user info from the cache or DB, or an error response.
    Does NOT add/update user in the DB.
    """
    try:
        # Check if the user information is already cached in the session
        if 'user' in session:
            return session['user'], None  # Return the user from session cache

        email = request.headers.get("X-Forwarded-User")
        persistent_id = request.headers.get("Shibboleth-Persistent-Id")

        if not email or not persistent_id:
            return None, (jsonify({"error": "Missing SSO authentication headers"}), 401)

        # Look up the user in your in-memory cache or DB by persistent_id
        user_data = user_cache.get(persistent_id)  # if using in-memory
        # or: user = User.query.filter_by(auth0_id=persistent_id).first()

        if not user_data:
            return None, (jsonify({"error": "SSO user not found"}), 403)

        # Cache the user data in session for future use
        session['user'] = user_data

        return user_data, None

    except Exception as e:
        logger.exception(f"validate_sso_user error: {e}")
        return None, (jsonify({"error": "Internal server error - user not validated"}), 500)


def maybe_require_user(ENABLE_AUTH):
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            if not ENABLE_AUTH:
                return f(user=None, *args, **kwargs)  # Inject `user=None` for consistency

            auth_method = session.get("auth_method")
            logger.info("------maybe_require_user----", auth_method)
            if not auth_method:
                return jsonify({"error": "Authentication method not set in session"}), 401

            if auth_method == "auth0":
                user, error_response = validate_and_get_user()
                if error_response:
                    logger.error(error_response)
                    return error_response
            elif auth_method == "shibboleth":
                user, error_response = validate_sso_user(request)
                if error_response:
                    logger.error(error_response)
                    return error_response
            else:
                return jsonify({"error": "Unknown authentication method"}), 400

            return f(user=user, *args, **kwargs)
        
        return decorated_function
    return decorator

