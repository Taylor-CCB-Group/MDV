import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def register_auth_routes(app):

    from flask import jsonify, request, redirect, url_for, session
    from mdvtools.dbutils.dbservice import UserService
    from mdvtools.dbutils.mdv_server_app import oauth
    from mdvtools.auth.authutils import update_cache
    from mdvtools.auth.auth0_provider import Auth0Provider
    

    """
    Registers the Auth0 routes like login, callback, logout, etc. to the Flask app,
    with centralized and route-specific error handling.
    """
    logger.info("Registering AUTH routes...")

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
                logger.info("auth0-login")
                session.clear()  
                return auth0_provider.login()
            except Exception as e:
                logger.exception(f"In register_auth0_routes : Error during login: {e}")
                return jsonify({"error": "Failed to start login process."}), 500

        # Route for the callback after login (handles the callback from Auth0)
        @app.route('/callback')
        def callback():
            try:
                logger.info("auth0-callback")
                code = request.args.get('code')  # Get the code from the callback URL
                if not code:
                    logger.error("Missing 'code' parameter in the callback URL.")
                    session.clear()  # Clear session if there's no code
                    return jsonify({"error": "Authorization code not provided."}), 400
                

                access_token = auth0_provider.handle_callback()
                if not access_token:  # If token retrieval fails, prevent redirecting
                    logger.error("Authentication failed: No valid token received.")
                    session.clear()  # Clear session in case of failure
                    return jsonify({"error": "Authentication failed."}), 401
                
            
                return redirect(url_for('index'))  # Redirect to the home page or any protected page
            except Exception as e:
                logger.exception(f"In register_auth0_routes : Error during callback: {e}")
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
                logger.exception(f"In register_auth0_routes: Error during logout: {e}")
                session.clear()
                return jsonify({"error": "Failed to log out."}), 500

        @app.route('/login_sso')
        def login_sso():
            """Redirect user to Shibboleth-protected login page on Apache."""
            try:
                # Extract SSO headers (available only *after* successful Shibboleth login)
                email = request.headers.get("X-Forwarded-User")
                persistent_id = request.headers.get("Shibboleth-Persistent-Id")
                
                # If headers not present, this is the *first* visit — redirect to Shibboleth login
                if not email or not persistent_id:
                    # Clear any existing session data to ensure we start with a fresh session
                    session.clear()
                    
                    # Store the authentication method as Shibboleth
                    session["auth_method"] = "shibboleth"  # Indicate Shibboleth login
                    session.modified = True
                    # Check if the Shibboleth login URL is provided in the environment
                    shibboleth_login_url = app.config.get('SHIBBOLETH_LOGIN_URL', None)

                    if shibboleth_login_url:
                        # Redirect the user to Shibboleth login page if the URL is configured
                        logger.info("Redirecting to Shibboleth login page...")
                        return redirect(shibboleth_login_url)
                    else:
                        # If Shibboleth URL is not provided, inform the user
                        logger.error("Shibboleth login URL not provided.")
                        return jsonify({"error": "Shibboleth login URL not provided."}), 500

                # User is authenticated — proceed to provision user
                user = UserService.add_or_update_user(
                    email=email,
                    auth0_id=persistent_id,
                    institution="University of Oxford"
                )

                user_data = {
                    "id": user.id,
                    "auth0_id": user.auth0_id,
                    "email": user.email,
                    "is_admin": user.is_admin
                }

                # Store user in cache
                update_cache(user_id=user.id, user_data=user_data)

                session["auth_method"] = "shibboleth"  # Indicate Shibboleth login
                session.modified = True

                logger.info(f"SSO login successful: {email}")
                return redirect("/")
            
            except Exception as e:
                # In case of error, clear the session and handle the error
                session.clear()  # Ensure session is cleared in case of failure
                logger.exception(f"In login_sso: Error during login: {e}")
                return jsonify({"error": "Failed to start login process using SSO."}), 500


        # You can also add a sample route to check the user's profile or token
        @app.route('/profile')
        def profile():
            try:
                auth_method = session.get('auth_method')

                if auth_method == 'auth0':
                    token = session.get('token')
                    if token:
                        user_info = auth0_provider.get_user(token)
                        if user_info:
                            return jsonify(user_info)
                        else:
                            return jsonify({"error": "Failed to retrieve user information from Auth0."}), 500
                    else:
                        return jsonify({"error": "Not authenticated with Auth0."}), 401

                elif auth_method == 'shibboleth':
                    # Retrieve user attributes from available Shibboleth headers
                    eppn = request.headers.get('Shibboleth-Eppn')  # e.g., whgu1064@ox.ac.uk
                    x_forwarded_user = request.headers.get('X-Forwarded-User')  # e.g., whgu1064@ox.ac.uk

                    if eppn:
                        user_data = {
                            "sub": eppn,
                            "first_name": "Unknown",
                            "last_name": "Unknown",
                            "email": eppn,
                            "association": "University of Oxford",
                            "avatarUrl": ""
                        }
                        return jsonify(user_data)
                    else:
                        return jsonify({"error": "Shibboleth attributes not found."}), 401

                else:
                    return jsonify({"error": "Unknown authentication method."}), 400

            except Exception as e:
                logger.exception(f"Error during profile retrieval: {e}")
                return jsonify({"error": "Failed to retrieve user profile."}), 500


        logger.info("Auth routes registered successfully!")

    except Exception as e:
        logger.exception(f"Error registering AUTH routes: {e}")
        raise