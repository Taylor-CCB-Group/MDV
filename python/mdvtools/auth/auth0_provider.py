import requests
from authlib.integrations.flask_client import OAuth
from flask import session, redirect, url_for
from typing import Optional
from auth_provider import AuthProvider

class Auth0Provider(AuthProvider):
    def __init__(self, app, oauth: OAuth, client_id: str, client_secret: str, domain: str):
        self.app = app
        self.oauth = oauth
        self.client_id = client_id
        self.client_secret = client_secret
        self.domain = domain
        self._initialize_oauth()

    def _initialize_oauth(self):
        self.oauth.register(
            'auth0',
            client_id=self.client_id,
            client_secret=self.client_secret,
            authorize_url=f'https://{self.domain}/authorize',
            authorize_params=None,
            access_token_url=f'https://{self.domain}/oauth/token',
            refresh_token_url=None,
            client_kwargs={'scope': 'openid profile email'},
        )

    def login(self) -> str:
        redirect_uri = url_for('callback', _external=True)
        return self.oauth.auth0.authorize_redirect(redirect_uri)

    def logout(self) -> None:
        session.clear()

    def get_user(self, token: str) -> Optional[dict]:
        user_info_url = f'https://{self.domain}/userinfo'
        response = requests.get(user_info_url, headers={'Authorization': f'Bearer {token}'})
        if response.status_code == 200:
            return response.json()
        return None

    def get_token(self) -> Optional[str]:
        return session.get('token')

    def handle_callback(self, code: str, redirect_uri: str) -> Optional[str]:
        token = self.oauth.auth0.authorize_access_token()
        session['token'] = token
        return token['access_token']

    def is_authenticated(self, token: str) -> bool:
        user_info = self.get_user(token)
        return user_info is not None
