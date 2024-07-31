""" import os
import json
from flask import Flask

# app = Flask(__name__, template_folder='../templates', static_folder='../static')
# static_folder='../../../dist/flask'
static_folder = "/app/dist/flask"
print(f">>>>> static path exists? {os.path.exists(static_folder)} <<<<<")
app = Flask(__name__, template_folder="../templates", static_folder=static_folder)

try:
    config_file_path = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "config.json"
    )
    print(config_file_path)
    with open(config_file_path) as config_file:
        config = json.load(config_file)
        app.config["SQLALCHEMY_DATABASE_URI"] = config.get("database_uri", "")
        app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = config.get(
            "track_modifications", False
        )
        app.config["upload_folder"] = config.get("upload_folder", "")
        app.config["projects_base_dir"] = config.get("projects_base_dir", "")
        print("In app.py : Configuration loaded successfully!")
except FileNotFoundError:
    print("Error: app.py script -> Configuration file not found.")
except Exception as e:
    print(f"An unexpected error occurred: {e}") """
