import os
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import ProjectBlueprint
from mdvtools.server import add_safe_headers
from flask import Flask, render_template, jsonify, request
import json
import threading
import random
import string

"""
Make a Flask app, and open all folders in ~/mdv/ as projects that can be served by it.
"""

# todo: make this configurable
# project_dir = os.path.join(os.path.expanduser("~"), "mdv")
project_dir = "/app/mdv"
# create the directory if it doesn't exist
if not os.path.exists(project_dir):
    os.makedirs(project_dir)

projects = [
    MDVProject(os.path.join(project_dir, d))
    for d in os.listdir(project_dir)
    if os.path.isdir(os.path.join(project_dir, d)) and d != 'lost+found'
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
    """watch the project folder for changes and update the projects list accordingly."""
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
            print(f"watcher adding '{p.id}'")
            try:
                p.serve(open_browser=False, app=app)
            except Exception as e:
                print(f"error serving {p.id}... {str(e)[:100]}")
    print("watcher exiting...")


if __name__ == "__main__":
    app = Flask(__name__)
    app.after_request(add_safe_headers)
    ProjectBlueprint.register_app(app)

    for p in projects:
        try:
            p.serve(open_browser=False, app=app)
        except Exception as e:
            print(f"error serving {p.id}... {str(e)[:100]}")

    @app.route("/")
    def index():
        # todo: figure out what to do here / how to configure routes
        return render_template("index.html")

    @app.route("/projects")
    def get_projects():
        # todo formalise relation between this, db-version of backend, frontend etc.
        return jsonify([{"id": p.id, "name": p.id} for p in projects])

    @app.route("/create_project", methods=["POST"])
    def create_project():
        try:
            project_id = (
                request.json["id"]
                if request.json and "id" in request.json
                else str("".join(random.choices(string.ascii_letters, k=6)))
            )
            # creating a project with the same id as an existing project is not allowed
            assert(project_id not in [p.id for p in projects])
            print(f"creating project '{project_id}'")
            p = MDVProject(os.path.join(project_dir, project_id))
            # in cases such as project_wizard, it's possible the project will already exist in this list
            # having been picked up by the watcher thread - so we should check for another project with the same id in `projects`
            if p.id not in [p.id for p in projects]:
                p.set_editable(True)
                projects.append(p)
                p.serve(app=app, open_browser=False)
            else:
                # the frontend doesn't care about this, and it's not an error given that 
                # we assert that the project doesn't already exist at the start of this function
                print(f"project '{p.id}' is already being served")
            return jsonify({"id": p.id, "name": p.id, "status": "success"})
        except Exception as e:
            return jsonify({"status": "error", "message": str(e)}), 500

    @app.route("/delete_project/<project_id>", methods=["DELETE"])
    def delete_project(project_id: str):
        """perhaps this should be routed by the project itself..."""
        try:
            # project_id = request.json["id"]
            print(f"deleting project '{project_id}'")
            p = next(p for p in projects if p.id == project_id)
            p.delete()
            projects.remove(p)
            return jsonify({"status": "success"})
        except Exception as e:
            return jsonify({"status": "error", "message": str(e)}), 500

    watcher = threading.Thread(target=watch_folder, args=(app,))
    # print("Oh frabjous day! Callooh! Callay!")
    watcher.daemon = True
    watcher.start()
    try:
        app.run(debug=True, port=5051)
    except Exception as e:
        print(e)
    running = False
