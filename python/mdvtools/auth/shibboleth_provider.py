from flask import session, redirect, request, jsonify
from mdvtools.dbutils.dbservice import UserService
from mdvtools.auth.authutils import update_cache
from mdvtools.auth.auth_provider import AuthProvider
import logging

logger = logging.getLogger(__name__)

class ShibbolethProvider(AuthProvider):
    def __init__(self, app):
        self.app = app
        self.logout_url = app.config.get("SHIBBOLETH_LOGOUT_URL")
        self.login_url = app.config.get("SHIBBOLETH_LOGIN_URL")

    def login(self):
        try:
            # Extract headers that Shibboleth sets after login
            email = request.headers.get("X-Forwarded-User")
            persistent_id = request.headers.get("Shibboleth-Persistent-Id")

            session["auth_method"] = "shibboleth"
            session.modified = True

            if not email or not persistent_id:
                session.clear()
                if self.login_url:
                    logger.info("Redirecting to Shibboleth login page...")
                    return redirect(self.login_url)
                else:
                    logger.error("Shibboleth login URL not configured.")
                    return jsonify({"error": "Shibboleth login URL not configured."}), 500

            # Provision user in DB
            user = UserService.add_or_update_user(
                email=email,
                auth_id=persistent_id,
                institution="University of Oxford"
            )

            # Cache user data
            user_data = {
                "id": user.id,
                "auth_id": user.auth_id,
                "email": user.email,
                "is_admin": user.is_admin
            }
            update_cache(user_id=user.id, user_data=user_data)

            logger.info(f"SSO login successful: {email}")
            return redirect("/")

        except Exception as e:
            session.clear()
            logger.exception(f"Error during Shibboleth login: {e}")
            return jsonify({"error": "Shibboleth login failed."}), 500

    def logout(self):
        session.clear()
        if self.logout_url:
            return redirect(self.logout_url)
        else:
            return jsonify({"error": "Shibboleth logout URL not configured."}), 500

    def get_user(self):
        eppn = request.headers.get('Shibboleth-Eppn')
        if eppn:
            return {
                "sub": eppn,
                "first_name": "Unknown",
                "last_name": "Unknown",
                "email": eppn,
                "association": "University of Oxford",
                "avatarUrl": ""
            }
        return None

    def get_token(self):
        return None  # Shibboleth uses headers

    def handle_callback(self):
        return None  # Not applicable
 
    def validate_user(self):
        """Validate user using Shibboleth."""
        from mdvtools.dbutils.dbmodels import User

        try:
            # Check if the user information is already cached in the session
            if 'user' in session:
                return session['user'], None  # Return the user from session cache

            email = request.headers.get("X-Forwarded-User")
            persistent_id = request.headers.get("Shibboleth-Persistent-Id")

            if not email or not persistent_id:
                return None, (jsonify({"error": "Missing SSO authentication headers"}), 401)

            # Look up user from database (or in-memory cache if applicable)
            user = User.query.filter_by(auth_id=persistent_id).first()
            if not user:
                return None, (jsonify({"error": "SSO user not found"}), 403)
            
            # Add the user to the in-memory cache
            user_data = {"id": user.id, "auth_id": user.auth_id, "email": user.email, "is_admin": user.is_admin}
            
            # Cache the user data in session for future use
            session['user'] = user_data
            session.modified = True

            return user_data, None

        except Exception as e:
            logger.exception(f"validate_sso_user error: {e}")
            return None, (jsonify({"error": "Internal server error - user not validated"}), 500)

    def sync_users_to_db(self) -> None:
        logger.info("Dummy sync - skipping user sync.")