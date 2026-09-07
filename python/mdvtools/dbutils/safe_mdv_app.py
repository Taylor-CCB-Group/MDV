from flask import Flask
import os
from mdvtools.logging_config import get_logger

# Setup logging
logger = get_logger(__name__)

logger.info("safe_mdv_app.py starting")

# https://stackoverflow.com/a/36033627/279703
class PrefixMiddleware(object):
    def __init__(self, app, prefix=''):
        self.app = app
        self.prefix = prefix

    def __call__(self, environ, start_response):
        if environ['PATH_INFO'].startswith(self.prefix):
            environ['PATH_INFO'] = environ['PATH_INFO'][len(self.prefix):]
            environ['SCRIPT_NAME'] = self.prefix
            return self.app(environ, start_response)
        else:
            response_body = b"This url does not belong to the app."
            headers = [('Content-Type', 'text/plain'), ('Content-Length', str(len(response_body)))]
            start_response('404 Not Found', headers)
            return [response_body]

# Import the app at the module level (outside any conditionals)
try:
    from mdvtools.dbutils.mdv_server_app import app
    if app is None:
        logger.error("app from mdv_server_app is None - using fallback")
        app = Flask(__name__)
except Exception as e:
    logger.exception(f'Error importing mdv_app: {e}')
    app = Flask(__name__)

# Apply prefix middleware at module level - works for both direct run and WSGI deployment
prefix = os.environ.get('MDV_API_ROOT', '')
if prefix:
    # Remove trailing slash from prefix if it exists
    if prefix.endswith('/'):
        prefix = prefix[:-1]
    
    logger.info(f"Applying PrefixMiddleware with prefix: {prefix}")
    app.wsgi_app = PrefixMiddleware(app.wsgi_app, prefix=prefix)

# Keep the existing __main__ block for direct execution
if __name__ == '__main__':
    """
    Running the app in local debugger with a path prefix
    This can be useful for testing behaviour app running on a sub-path
    """
    from mdvtools.websocket import socketio
    from mdvtools.socketio_upload import initialize_socketio_upload
    
    logger.info("running as __main__")
    # Set default for local testing if not already set
    if 'MDV_API_ROOT' not in os.environ:
        os.environ['MDV_API_ROOT'] = '/test/'

    # Register upload socket events
    if initialize_socketio_upload is not None:
        initialize_socketio_upload(app)
    
    # Run with socketio if available, otherwise fallback to Flask
    if socketio is not None:
        socketio.run(app, host='0.0.0.0', debug=False, port=5056)
    else:
        app.run(host='0.0.0.0', debug=False, port=5056)