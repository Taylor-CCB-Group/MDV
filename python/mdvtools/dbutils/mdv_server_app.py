import os
import json
# import threading
# import random
# import string
# from flask import Flask, render_template, jsonify, request
from mdvtools.server import add_safe_headers
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprint
from mdvtools.dbutils.app import app
from mdvtools.dbutils.dbmodels import db, Project
from mdvtools.dbutils.routes import register_global_routes

def create_projects_from_filesystem(base_dir):
    try:
        print(f"Scanning base directory: {base_dir}")

        # Get all project paths from the database
        projects_in_db = {project.path for project in Project.query.with_entities(Project.path).all()}
        print(f"Project paths in DB: {projects_in_db}")

        # Get all project directories in the filesystem
        project_paths_in_fs = {os.path.join(base_dir, d) for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))}
        print(f"Project paths in filesystem: {project_paths_in_fs}")

        # Determine which project paths are in the filesystem but not in the database
        missing_project_paths = project_paths_in_fs - projects_in_db
        print(f"Missing project paths: {missing_project_paths}")

        # Iterate over missing project paths to create and serve them
        for project_path in missing_project_paths:
            print(f"Processing project path: {project_path}")
            
            if os.path.exists(project_path):
                try:
                    project_name = os.path.basename(project_path)

                    p = MDVProject(project_path)
                    p.set_editable(True)
                    p.serve(app=app, open_browser=False)
                    print(f"Serving project: {project_path}")

                    # Create a new Project record in the database with the default name
                    new_project = Project(name=project_name, path=project_path)
                    db.session.add(new_project)
                    db.session.commit()
                    print(f"Added project to DB: {new_project}")
                except Exception as e:
                    print(f"In create_projects_from_filesystem: Error creating project at path '{project_path}': {e}")
            else:
                print(f"In create_projects_from_filesystem: Error - Project path '{project_path}' does not exist.")
                
    except Exception as e:
        print(f"In create_projects_from_filesystem: Error retrieving projects from database: {e}")

if __name__ == '__main__':
    try:
        config_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.json')
        with open(config_file_path) as config_file:
            config = json.load(config_file)
            base_dir = config.get('projects_base_dir', 'mdv')
            print(f"Loaded config: {config}")
    except FileNotFoundError:
        print("Error: mdv_server_app.py -> Configuration file not found.")
        exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while loading configuration: {e}")
        exit(1)

    if not os.path.exists(base_dir):
        try:
            os.makedirs(base_dir)
            print(f"Created base directory: {base_dir}")
        except Exception as e:
            print(f'Error creating base directory: {e}')
            exit(1)

    app.after_request(add_safe_headers)

    with app.app_context():
        try:
            print("Initialized app context")
            
            db.create_all()
            print("Created the database tables")

            print("Registering the global routes")
            register_global_routes(base_dir)
            
            print("Registering the blueprint(register_app)")
            ProjectBlueprint.register_app(app)
            
            print("Start- create_projects_from_filesystem")
            create_projects_from_filesystem(base_dir)

        except Exception as e:
            print(f'Error initializing app: {e}')

    try:
        app.run(host='0.0.0.0', debug=True, port=5055)
    except Exception as e:
        print(f'Error running app: {e}')
