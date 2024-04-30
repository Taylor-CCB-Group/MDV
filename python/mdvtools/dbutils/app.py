import os
import json
from flask import Flask

app = Flask(__name__)

try:
    config_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.json')
    print(config_file_path)
    with open(config_file_path) as config_file:
        config = json.load(config_file)
        app.config['SQLALCHEMY_DATABASE_URI'] = config.get('database_uri', '')
        app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = config.get('track_modifications', False)
        # app.config['UPLOAD_FOLDER'] = config.get('upload_folder', '')
        print("Configuration loaded successfully!")
except FileNotFoundError:
    print("Error: app.py script -> Configuration file not found.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
