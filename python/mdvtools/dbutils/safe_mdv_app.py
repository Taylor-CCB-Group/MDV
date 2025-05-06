from flask import Flask
import logging

# Setup logging
logger = logging.getLogger(__name__)

logger.info("safe_mdv_app.py starting")

# Import the app at the module level (outside any conditionals)
try:
    from mdvtools.dbutils.mdv_server_app import app
    logger.info("imported mdv_app at module level")
except Exception as e:
    logger.exception(f'Error importing mdv_app: {e}')
    app = Flask(__name__)

# Keep the existing __main__ block for direct execution
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logger.info("running as __main__")
    app.run(host='0.0.0.0', debug=False, port=5055)