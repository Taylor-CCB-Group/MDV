import os

# mdv_server_app builds the app and scans the project root when it is imported.
# These tests import functions from it and must not touch the real deployment.
os.environ["MDV_SKIP_SERVER_STARTUP"] = "1"
