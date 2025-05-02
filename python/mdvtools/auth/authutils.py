from flask import jsonify, session, request, current_app
from functools import wraps

def validate_and_get_user():
    from mdvtools.auth.auth0_provider import Auth0Provider
    from mdvtools.dbutils.mdv_server_app import oauth
    from mdvtools.dbutils.mdv_server_app import user_cache
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
            return None, jsonify({"error": "Authentication required"}), 401

        # Validate token
        if not auth0_provider.is_token_valid(token):
            return None, jsonify({"error": "Invalid or expired token"}), 401

        # Retrieve user info from Auth0
        user_info = auth0_provider.get_user({"access_token": token})

        if not user_info:
            return None, jsonify({"error": "User not found"}), 404

        # Get Auth0 user ID
        auth0_id = user_info.get("sub")

        # Fetch user directly from in-memory cache
        cached_user = user_cache.get(auth0_id)
        if cached_user:
            session['user'] = cached_user
            return cached_user, None  # Return cached user

        # This should rarely happen since all users were cached at app startup.
        # But if the user is missing from Redis, fallback to DB and log the issue.
        print(f"User {auth0_id} not found in cache, falling back to DB!")

        # Query the user from the database if not in cache
        user = User.query.filter_by(auth0_id=auth0_id).first()
        if not user:
            return None, jsonify({"error": "User not found"}), 404

        # Add the user to the in-memory cache
        user_data = {"id": user.id, "auth0_id": user.auth0_id, "email": user.email, "is_admin": user.is_admin}
        user_cache[auth0_id] = user_data

        # Cache the user data in session for future use
        session['user'] = user_data

        return user_data, None

    except Exception as e:
        print(f"Error during authentication: {e}")
        return None, jsonify({"error": "Internal server error - user not validated"}), 500
    

def validate_sso_user(request):
    from mdvtools.dbutils.mdv_server_app import user_cache
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
            return None, jsonify({"error": "Missing SSO authentication headers"}), 401

        # Look up the user in your in-memory cache or DB by persistent_id
        user_data = user_cache.get(persistent_id)  # if using in-memory
        # or: user = User.query.filter_by(auth0_id=persistent_id).first()

        if not user_data:
            return None, jsonify({"error": "SSO user not found"}), 403

        # Cache the user data in session for future use
        session['user'] = user_data

        return user_data, None

    except Exception as e:
        print(f"validate_sso_user error: {e}")
        return None, jsonify({"error": "Internal server error - user not validated"}), 500


def maybe_require_user(ENABLE_AUTH):
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            if not ENABLE_AUTH:
                return f(user=None, *args, **kwargs)  # Inject `user=None` for consistency

            auth_method = session.get("auth_method")
            print("------maybe_require_user----", auth_method)
            if not auth_method:
                return jsonify({"error": "Authentication method not set in session"}), 401

            if auth_method == "auth0":
                user, error_response = validate_and_get_user()
                if error_response:
                    print("------maybe_require_user----ERROR in auth0 validation")
                    return error_response
            elif auth_method == "shibboleth":
                user, error_response = validate_sso_user(request)
                if error_response:
                    print("------maybe_require_user----ERROR in sso validation")
                    return error_response
            else:
                return jsonify({"error": "Unknown authentication method"}), 400

            return f(user=user, *args, **kwargs)
        
        return decorated_function
    return decorator

