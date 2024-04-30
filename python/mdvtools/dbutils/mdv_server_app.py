import os
import json
from mdvtools.server import add_safe_headers
from mdvtools.dbutils.app import app
from mdvtools.dbutils.dbmodels import db
from mdvtools.dbutils.routes import register_routes
from mdvtools.dbutils.project_loader import create_all_projects


def serve_projects(projects):
    try:
        register_routes(projects)

        for p in projects:
            try:
                print(p.name)
                p.serve(open_browser=False, app=app)
            except Exception as e:
                print(f'Error serving {p.name}: {e}')
    except Exception as e:
        print(f'Error registering routes: {e}')


if __name__ == '__main__':

    try:
        config_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.json')
        with open(config_file_path) as config_file:
            config = json.load(config_file)
            base_dir = config.get('projects_base_dir', 'mdv')
    except FileNotFoundError:
        print("Error: mdv_server_app.py -> Configuration file not found.")
        exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while loading configuration: {e}")
        exit(1)

    #base_dir = os.path.join(os.path.expanduser('~'), 'mdv')

    if not os.path.exists(base_dir):
        try:
            os.makedirs(base_dir)
        except Exception as e:
            print(f'Error creating base directory: {e}')
            exit(1)

    app.after_request(add_safe_headers)

    with app.app_context():
        try:
            db.create_all()
            projects = create_all_projects(base_dir)
            serve_projects(projects)
        except Exception as e:
            print(f'Error initializing app: {e}')

    try:
        app.run(debug=True, port=5051)
    except Exception as e:
        print(f'Error running app: {e}')
