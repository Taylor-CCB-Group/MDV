from flask import (
    Flask,
    Blueprint,
    render_template,
    request,
    make_response,
    send_file,
    Response,
    jsonify,
)
import webbrowser
import mimetypes
import json
import sys
import re
from werkzeug.security import safe_join
from mdvtools.websocket import mdv_socketio
from mdvtools.mdvproject import MDVProject
import os
import pandas as pd
from typing import Optional
from datetime import datetime

routes = set()


# consider using flask_cors...
def add_safe_headers(resp):
    # headers required for web workers
    resp.headers["Cross-Origin-Opener-Policy"] = "same-origin"
    resp.headers["Cross-Origin-Embedder-Policy"] = "require-corp"
    # headers required if serving endpoints for another server e,g dev server
    resp.headers["Access-Control-Allow-Origin"] = "*"
    resp.headers["Access-Control-Allow-Headers"] = "Content-Type"
    return resp


# flask send_file can't always cope with relative paths
# sets the cwd to the python path for some reason
def _send_file(f):
    if not os.path.isabs(f):
        f = os.path.join(os.getcwd(), f)
    return send_file(f)


def get_range(file_name, range_header):
    file = open(file_name, "rb")
    size = sys.getsizeof(file_name)
    byte1, byte2 = 0, None

    m = re.search(r"(\d+)-(\d*)", range_header)
    if not m:
        raise Exception("Invalid Range Header")
    g = m.groups()

    if g[0]:
        byte1 = int(g[0])
    if g[1]:
        byte2 = int(g[1])

    length = size - byte1
    if byte2 is not None:
        length = byte2 - byte1 + 1

    file.seek(byte1)
    data = file.read(length)
    rv = Response(
        data, 206, mimetype=mimetypes.guess_type(file_name)[0], direct_passthrough=True
    )
    rv.headers.add(
        "Content-Range", "bytes {0}-{1}/{2}".format(byte1, byte1 + length - 1, size)
    )
    rv.headers.add("Accept-Ranges", "bytes")
    file.close()
    return rv


def create_app(
    project: MDVProject,
    open_browser=True,
    port=5050,
    websocket=False,
    app: Optional[Flask] = None,
):
    if app is None:
        route = ""
        # route = "/project/" + project.name # for testing new API with simple app...
        app = Flask(__name__)
        print(f"created Flask {app}")
        # add headers to allow web workers
        app.after_request(add_safe_headers)
        project_bp = app
        multi_project = False
        # nb, may make this default to False again soon.
        ### 'MEW' in Unity is using IWeb PostMessage, not WebSockets.
        ### but this will be used for local testing in short-term, and potentially other things later.
        if websocket:
            mdv_socketio(app)
    else:
        # add routes for this project to existing app
        # set the route prefix to the project name, derived from the dir name.
        # this is to allow multiple projects to be served from the same server.
        multi_project = True
        route = "/project/" + project.name + "/"
        project_bp = Blueprint(project.name, __name__, url_prefix=route)
    if route in routes:
        raise Exception(
            "Route already exists - can't have two projects with the same name"
        )
    routes.add(route)

    @project_bp.route("/")
    def project_index():
        print("recieved request to project_index")
        # `_mdvInit('{{route}}')` in template...
        return render_template("page.html", route=route)

    @project_bp.route("/<file>.b")
    def get_binary_file(file):
        # should this now b '.gz'?
        file_name = safe_join(project.dir, file + ".b")
        range_header = request.headers.get("Range", None)
        return get_range(file_name, range_header)

    # duplicate of above, but for .gz files in case that's needed.
    # (there was some reason for changing to this, but I can't fully remember the status
    # so maybe better to support both for now)
    @project_bp.route("/<file>.gz")
    def get_binary_file_gz(file):
        file_name = safe_join(project.dir, file + ".gz")
        range_header = request.headers.get("Range", None)
        return get_range(file_name, range_header)

    @project_bp.route("/<file>.json")
    def get_json_file(file: str):
        if project.dir is None:
            return "Project directory not found", 404
        path = safe_join(project.dir, file + ".json")
        if path is None or not os.path.exists(path):
            return "File not found", 404
        return send_file(path)

    # empty page to put popout content
    @project_bp.route("/popout.html")
    def popout():
        return "<!DOCTYPE html><html><head></head><body></body></html>"

    # gets the raw byte data and packages it in the correct response
    @project_bp.route("/get_data", methods=["POST"])
    def get_data():
        try:
            data = request.json
            if not data or "columns" not in data or "data_source" not in data:
                raise Exception(
                    "Request must contain JSON with 'columns' and 'data_source'"
                )
            bytes_ = project.get_byte_data(data["columns"], data["data_source"])
            response = make_response(bytes_)
            response.headers.set("Content-Type", "application/octet-stream")
            return response
        except Exception as e:
            print(e)
            return "Problem handling request", 400

    # images contained in the project
    @project_bp.route("/images/<path:path>")
    def images(path):
        try:
            return _send_file(project.get_image(path))
        except Exception:
            return _send_file(safe_join(project.imagefolder, path))

    # All the project's metadata
    @project_bp.route("/get_configs", methods=["GET", "POST"])
    def get_configs():
        return jsonify(project.get_configs())

    # gets a particular view
    @project_bp.route("/get_view", methods=["POST"])
    def get_view():
        data = request.json
        if not data or "view" not in data:
            return "Request must contain JSON with 'view'", 400
        return jsonify(project.get_view(data["view"]))

    # get any custom row data
    @project_bp.route("/get_row_data", methods=["POST"])
    def get_row_data():
        req = request.json
        if req is None:
            return json.dumps({"data": None})
        path = safe_join(
            project.dir, "rowdata", req["datasource"], f"{req['index']}.json"
        )
        if path is None or not os.path.exists(path):
            return json.dumps({"data": None})
        with open(path) as f:
            if f is None:
                return json.dumps({"data": None})
            return f.read()

    # get arbitrary data
    @project_bp.route("/get_binary_data", methods=["POST"])
    def get_binary_data():
        req = request.json
        try:
            if req is None or "datasource" not in req or "name" not in req:
                return "Request must contain JSON with 'datasource' and 'name'", 400
            if project.dir is None or not os.path.exists(project.dir):
                return "Project directory not found", 404
            path = safe_join(
                project.dir, "binarydata", req["datasource"], f"{req['name']}.gz"
            )
            if path is None or not os.path.exists(path):
                return "Binary data not found", 404
            with open(path, "rb") as f:
                data = f.read()
        except Exception:
            # data='' # satisfy type checker - was None, haven't tested if this is better or worse.
            # probably better to return an error.
            return "Problem getting binary data", 500
        return data

    # only the specified region of track files (bam,bigbed,tabix)
    # needs to be returned
    @project_bp.route("/tracks/<path:path>")
    def send_track(path):
        file_name = safe_join(project.trackfolder, path)
        range_header = request.headers.get("Range", None)
        if not range_header:
            return _send_file(file_name)
        return get_range(file_name, range_header)

    @project_bp.route("/save_state", methods=["POST"])
    def save_data():
        success = True
        try:
            state = request.json
            project.save_state(state)
        except Exception:
            success = False

        return jsonify({"success": success})

    @project_bp.route("/add_datasource", methods=["POST"])
    def add_datasource():
        if (
            "permission" not in project.state
            or not project.state["permission"] == "edit"
        ):
            return "Project is read-only", 400
        success = True
        try:
            name = request.form["name"]
            if not name:
                return "Request must contain 'name'", 400
            # xxx - not how column metadata should be passed, todo fix
            # cols = (
            #     request.form["columns"].split(",")
            #     if "columns" in request.form
            #     else None
            # )
            view = request.form["view"] if "view" in request.form else None
            replace = True if "replace" in request.form else False
            if not replace and name in [ds['name'] for ds in project.datasources]:
                return f"Datasource '{name}' already exists, and 'replace' was not set in request", 400
            if "file" not in request.files:
                return "No 'file' provided in request form data", 400
            file = request.files["file"]
            supplied_only = True if "supplied_only" in request.form else False
            if not file or file.mimetype != "text/csv":
                return "File must be a CSV", 400
            file.seek(0)
            # will this work? can we return progress to the client?
            df = pd.read_csv(file.stream)
            project.add_datasource(
                name,
                df,
                # cols,
                add_to_view=view,
                supplied_columns_only=supplied_only,
                replace_data=replace,
            )
            if "backend" in request.form:
                response = add_datasource_backend(project,view)
                return response
        except Exception as e:
            # success = False
            return str(e), 400

        metadata = project.get_datasource_metadata(name)
        return jsonify({"success": success, "metadata": metadata})

    if open_browser:
        webbrowser.open(f"http://localhost:{port}/{route}")

    if multi_project:
        if not isinstance(project_bp, Blueprint):
            raise Exception(
                "assert: project_bp must be a Flask instance when multi_project is True"
            )
        print(f"Adding project {project.name} to existing app")
        app.register_blueprint(project_bp)
    else:
        app.run(host="0.0.0.0", port=port, debug=True)

def add_datasource_backend(project,view):
    from mdvtools.dbutils.dbmodels import db, Project, File, User
    try:
        project_name = project.name
        if not project_name:
            return jsonify({'error': 'Project name is missing.'}), 400
        
        # Check if project exists
        project_db = Project.query.filter_by(name=project_name).first()
        if not project_db:
            return jsonify({'project does not exist in database'}), 400

        
        file_set = [project.h5file, project.datasourcesfile]
        if view:
                file_set.append(project.viewsfile)

        for file in file_set:
            existing_file = File.query.filter_by(name=os.path.basename(file), project_id=project_db.id).first()
            if existing_file:
               
                # Update the database entry with new file path and update timestamp
               
                existing_file.update_timestamp = datetime.now()

                

            else:
                
                # Create a new file entry
                new_file = File(name=os.path.basename(file), file_path=None, project=project_db) # type: ignore
                db.session.add(new_file)
                
        db.session.commit()
        return jsonify({'message': f'Files {file_set} have been modified under project "{project_name}  "'}), 200
    except Exception as e:
        return jsonify({'error': str(e)}), 500