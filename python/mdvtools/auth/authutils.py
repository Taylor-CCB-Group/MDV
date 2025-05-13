import logging
from flask import session, request,redirect, current_app

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# in_memory_cache.py
user_cache = {}  # key: auth_id -> user details
user_project_cache = {}  # key: user_id -> project permissions
all_users_cache = []  # list of all user summaries
active_projects_cache = []


def get_auth_provider():
    # Try to fetch the authentication method from session first
    auth_method = session.get("auth_method", None)
    
    if not auth_method:
        # If no auth method is found in the session, use the default from app config
        auth_method = current_app.config.get("DEFAULT_AUTH_METHOD", "auth0")
    
    if not auth_method:
        raise ValueError("Authentication method is not specified in session or config.")
    
    
    auth_method = auth_method.lower()

    if auth_method == "auth0":
        from mdvtools.auth.auth0_provider import Auth0Provider
        return Auth0Provider(
            current_app,
            oauth=current_app.extensions.get('oauth'),
            client_id=current_app.config['AUTH0_CLIENT_ID'],
            client_secret=current_app.config['AUTH0_CLIENT_SECRET'],
            domain=current_app.config['AUTH0_DOMAIN']
        )

    elif auth_method == "shibboleth":
        from mdvtools.auth.shibboleth_provider import ShibbolethProvider
        return ShibbolethProvider(current_app)

    elif auth_method == "dummy":
        from mdvtools.auth.dummy_provider import DummyAuthProvider
        return DummyAuthProvider(current_app)
    else:
        raise ValueError(f"Unsupported auth method: {auth_method}")

def is_authenticated():
    """Validate the current user via Auth0 or Shibboleth."""
    ENABLE_AUTH = current_app.config.get("ENABLE_AUTH", False)
    
    if not ENABLE_AUTH:
        return True
    
    try:
        provider = get_auth_provider()
    except Exception as e:
        current_app.logger.error(f"Failed to get auth provider: {e}")
        return False

    user, error_response = provider.validate_user()
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


def cache_user_projects():
    """
    Caches user details and their associated project permissions in memory.
    """

    from mdvtools.dbutils.dbmodels import User, UserProject
    from mdvtools.dbutils.dbservice import ProjectService
    
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
                "auth_id": user.auth_id,
                "email": user.email,
                "is_admin": user.is_admin
            }
            user_cache[user.auth_id] = user_data
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

        logger.info("Cached users and their project permissions in memory.")
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



