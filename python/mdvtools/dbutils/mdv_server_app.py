import os
import json
#import threading
#import random
#import string
#from flask import Flask, render_template, jsonify, request
from mdvtools.server import add_safe_headers
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprint
from mdvtools.dbutils.app import app
from mdvtools.dbutils.dbmodels import db, Project
from mdvtools.dbutils.routes import register_global_routes



def create_projects_from_filesystem(base_dir):

    try:
        # Get all project IDs from the database
        project_ids_in_db = {project.id for project in Project.query.with_entities(Project.id).all()}

        # Get all project directories in the filesystem (assumed to be digits only)
        project_ids_in_fs = {int(d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)) and d.isdigit()}

        missing_project_ids = project_ids_in_fs - project_ids_in_db

        for project_id in missing_project_ids:
            project_path = os.path.join(base_dir, str(project_id))
            if os.path.exists(project_path):
                try:
                    p = MDVProject(project_path)
                    p.set_editable(True)
                    p.serve(app=app, open_browser=False)
                    
                    # Create a new Project record in the database with the default name
                    new_project = Project()
                    db.session.add(new_project)
                    db.session.commit()
                except Exception as e:
                    print(f"In create_projects_from_filesystem: Error creating project with ID '{project_id}': {e}")
            else:
                print(f"In create_projects_from_filesystem:  Error- Project path '{project_path}' does not exist.")
                

    except Exception as e:
        print(f"In create_projects_from_filesystem: Error retrieving projects from database: {e}")



# running = True
# def watch_folder(app: Flask):
#     """watch the project folder for changes and update the projects list accordingly."""
#     import time

#     while running:
#         time.sleep(2)
#         existing_project_dirs = [p.dir for p in projects]
#         new_projects = [
#             MDVProject(os.path.join(project_dir, d))
#             for d in os.listdir(project_dir)
#             if os.path.isdir(os.path.join(project_dir, d))
#             and os.path.join(project_dir, d) not in existing_project_dirs
#         ]
#         projects.extend(new_projects)
#         for p in new_projects:
#             print(f"watcher adding '{p.id}'")
#             try:
#                 p.serve(open_browser=False, app=app)
#             except Exception:
#                 print(f"error serving {p.id}...")
#     print("watcher exiting...")


if __name__ == "__main__":
    try:
        config_file_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "config.json"
        )
        with open(config_file_path) as config_file:
            config = json.load(config_file)
            base_dir = config.get("projects_base_dir", "mdv")
    except FileNotFoundError:
        print("Error: mdv_server_app.py -> Configuration file not found.")
        exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while loading configuration: {e}")
        exit(1)

    # base_dir = os.path.join(os.path.expanduser('~'), 'mdv')

    if not os.path.exists(base_dir):
        try:
            os.makedirs(base_dir)
        except Exception as e:
            print(f"Error creating base directory: {e}")
            exit(1)

    app.after_request(add_safe_headers)

    with app.app_context():
        try:
            db.create_all()
            register_global_routes(base_dir)
            ProjectBlueprint.register_app(app)
            
            create_projects_from_filesystem(base_dir)

            
            
        except Exception as e:
            print(f'Error initializing app: {e}')

    # watcher = threading.Thread(target=watch_folder, args=(app,))
    # # print("Oh frabjous day! Callooh! Callay!")
    # watcher.daemon = True
    # watcher.start()
    
    try:
        app.run(host="0.0.0.0", debug=True, port=5055)
    except Exception as e:
        print(f"Error running app: {e}")
