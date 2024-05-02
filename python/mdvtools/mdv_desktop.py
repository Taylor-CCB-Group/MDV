import os
from mdvtools.mdvproject import MDVProject
from .server import add_safe_headers
from flask import Flask, render_template, jsonify
import json
import threading

"""
Make a Flask app, and open all folders in ~/mdv/ as projects that can be served by it.
"""

# todo: make this configurable
project_dir = os.path.join(os.path.expanduser("~"), "mdv")
# create the directory if it doesn't exist
if not os.path.exists(project_dir):
    os.makedirs(project_dir)

projects = [
    MDVProject(os.path.join(project_dir, d))
    for d in os.listdir(project_dir)
    if os.path.isdir(os.path.join(project_dir, d))
]
# add other projects as listed in config file (maybe later a sqlite db or something)
# nothing more permanent than the provisional... just want a quick way of including projects on an external volume.
# manually editted for now, and not watched for changes.
config_file = os.path.join(project_dir, "config.json")
if os.path.exists(config_file):
    import json

    with open(config_file, "r") as f:
        config = json.load(f)
    for d in config["projects"]:
        print(f"adding project '{d}' from config file")
        try:
            projects.append(MDVProject(d))
        except Exception:
            print(f"error adding project '{d}' from config file")
else:
    with open(config_file, "w") as f:
        json.dump({"projects": []}, f)

running = True
def watch_folder(app: Flask):
    '''watch the project folder for changes and update the projects list accordingly.'''
    import time

    while running:
        time.sleep(2)
        existing_project_dirs = [p.dir for p in projects]
        new_projects = [
            MDVProject(os.path.join(project_dir, d))
            for d in os.listdir(project_dir)
            if os.path.isdir(os.path.join(project_dir, d))
            and os.path.join(project_dir, d) not in existing_project_dirs
        ]
        projects.extend(new_projects)
        for p in new_projects:
            print(f"watcher adding '{p.name}'")
            try:
                p.serve(open_browser=False, app=app)
            except Exception:
                print(f"error serving {p.name}...")
    print("watcher exiting...")

if __name__ == "__main__":
    app = Flask(__name__)
    app.after_request(add_safe_headers)

    for p in projects:
        try:
            p.serve(open_browser=False, app=app)
        except Exception:
            print(f"error serving {p.name}...")

    @app.route("/")
    def index():
        # todo: figure out what to do here / how to configure routes
        return render_template("index.html")

    @app.route("/projects")
    def get_projects():
        return jsonify([p.name for p in projects])

    watcher = threading.Thread(target=watch_folder, args=(app,))
    watcher.start()
    app.run(debug=True, port=5051)
    running = False
