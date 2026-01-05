import os
from mdvtools.logging_config import get_logger

logger = get_logger(__name__)

def register_auth_routes(app):

    from flask import jsonify, request, redirect, url_for, session
    from mdvtools.auth.authutils import get_auth_provider
    
    """
     Registers the Auth0 routes like login, callback, logout, etc. to the Flask app,
    with centralized and route-specific error handling.
    """
    logger.info("Registering AUTH routes...")

    try:
        # Route for login (redirects to Auth0 for authentication)
        @app.route('/login')
        def login():
            """Handles standard login (e.g., Auth0)."""
            try:
                logger.info("Login initiated")
                session.clear()
                session["auth_method"] = app.config.get("DEFAULT_AUTH_METHOD", "dummy").lower()
                session.modified = True
                return get_auth_provider().login()
            except Exception as e:
                session.clear()
                logger.exception(f"In register_auth_routes : Error during login: {e}")
                return jsonify({"error": "Failed to start login process."}), 500

        # Route for the callback after login (handles the callback from Auth0)
        @app.route('/callback')
        def callback():
            try:
                logger.info("Callback route hit")
                code = request.args.get('code')

                # If this is a code-based provider like Auth0, validate presence of 'code'
                if not code and session.get("auth_method", "").lower() == "auth0":
                    logger.error("Missing 'code' parameter in the callback URL for Auth0.")
                    session.clear()
                    return jsonify({"error": "Authorization code not provided."}), 400

                auth_provider = get_auth_provider()

                access_token = auth_provider.handle_callback()
                if not access_token:
                    logger.error("Authentication failed: No valid token received.")
                    session.clear()
                    return jsonify({"error": "Authentication failed."}), 401

                logger.info("Authentication successful. Redirecting based on configuration.")
                # Check for MDV_API_ROOT environment variable
                mdv_api_root = os.getenv('MDV_API_ROOT')
                if mdv_api_root:
                    logger.info(f"Redirecting to MDV_API_ROOT: {mdv_api_root}")
                    return redirect(mdv_api_root)
                else:
                    logger.info("MDV_API_ROOT not set, redirecting to index")
                    return redirect(url_for("index"))
            except Exception as e:
                logger.exception(f"In register_auth_routes : Error during callback: {e}")
                session.clear()  # Clear session on error
                return jsonify({"error": "Failed to complete authentication process."}), 500

        # Route for logout (clears the session and redirects to home)
        @app.route('/logout')
        def logout():
            try:
                auth_method = session.get("auth_method")
                logger.info(f"Logout initiated for auth method: {auth_method}")

                # Always attempt to clear the session early
                session.clear()
                session.modified = True

                # Get the appropriate provider if method is set
                if auth_method:
                    auth_provider = get_auth_provider()
                    return auth_provider.logout()

                # Fallback for unknown or missing auth method
                logger.warning("No auth method found in session during logout.")
                return jsonify({"message": "Session cleared, but auth method unknown."}), 200

            except Exception as e:
                logger.exception(f"In register_auth_routes: Error during logout: {e}")
                session.clear()
                return jsonify({"error": "Failed to log out."}), 500

        @app.route('/login_sso')
        def login_sso():
            """Redirect user to Shibboleth-protected login page on Apache."""
            try:
                logger.info("Shibboleth SSO login initiated")
                session.clear()  # Clear previous session
                session["auth_method"] = "shibboleth"
                session.modified = True
                return get_auth_provider().login()
            except Exception as e:
                # In case of error, clear the session and handle the error
                session.clear()  # Ensure session is cleared in case of failure
                logger.exception(f"In register_auth_routes: Error during login SSO: {e}")
                return jsonify({"error": "Failed to start login process using SSO."}), 500


        # You can also add a sample route to check the user's profile or token
        @app.route('/profile')
        def profile():
            try:
                # Get the appropriate authentication provider based on session
                provider = get_auth_provider()
                # If the provider is Auth0, we need to pass the token (which should be in the session)
                token = session.get("token")

                # Get the user information
                if token:
                    # If token is available, use it for Auth0
                    user_info = provider.get_user(token=token)
                else:
                    # If no token (e.g., using Shibboleth), just call get_user without token
                    user_info = provider.get_user()

                if not user_info:
                    return jsonify({"error": "User profile could not be retrieved."}), 404
                # Return the user profile in a JSON response
                return jsonify(user_info), 200
            except Exception as e:
                logger.exception(f"Error during profile retrieval: {e}")
                return jsonify({"error": "Failed to retrieve user profile."}), 500


        logger.info("Auth routes registered successfully!")

    except Exception as e:
        logger.exception(f"Error registering AUTH routes: {e}")
        raise