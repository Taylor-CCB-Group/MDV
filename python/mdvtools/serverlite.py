"""
Lightweight Flask server for serving MDV project data and metadata.

This module provides a minimal Flask application that exposes endpoints
for serving project data, images, configurations, and handling file operations
with support for HTTP range requests and security headers.
"""
import json
import os
import sys
from flask import (
    Flask,
    render_template,
    request,
    make_response,
    jsonify,
)

from werkzeug.security import safe_join
from mdvtools.mdvproject import MDVProject

from mdvtools.server_utils import (
    send_file,
    get_range,
    add_safe_headers,
)

def create_app(project: MDVProject):
    """
    Create and configure a Flask application for serving MDV project data.

    Args:
        project (MDVProject): The MDV project instance to serve

    Returns:
        Flask: Configured Flask application
    """
    app = Flask(__name__)
    app.after_request(add_safe_headers)

    @app.route("/")
    def project_index():
        return render_template("page.html")

    # gets the raw byte data and packages it in the correct response
    @app.route("/get_data", methods=["POST"])
    def get_data():
        try:
            data = request.json
            if not data or "columns" not in data or "data_source" not in data:
                raise ValueError(
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
    @app.route("/images/<path:path>")
    def images(path):
        response = make_response(send_file(safe_join(project.imagefolder, path)))
        response.headers['Cache-Control'] = 'public, max-age=31536000'  # Cache for 1 year
        return response

    # All the project's metadata
    # this method is no longer used - configs are referenced quite directly
    @app.route("/get_configs", methods=["GET", "POST"])
    def get_configs():
        return jsonify(project.get_configs())
    
    # json files are obtained by url instead  of 
    @app.route("/<file>.json")
    def get_json_file(file: str):
        if project.dir is None:
            return "Project directory not found", 404
        path = safe_join(project.dir, file + ".json")
        if path is None or not os.path.exists(path):
            return "File not found", 404
        return send_file(path)


    # gets a particular view
    @app.route("/get_view", methods=["POST"])
    def get_view():
        data = request.json
        if not data or "view" not in data:
            return "Request must contain JSON with 'view'", 400
        return jsonify(project.get_view(data["view"]))

    # get any custom row data
    @app.route("/get_row_data", methods=["POST"])
    def get_row_data():
        req = request.json
        if req is None:
            return json.dumps({"data": None})
        path = safe_join(
            project.dir, "rowdata", req["datasource"], f"{req['index']}.json"
        )
        if path is None or not os.path.exists(path):
            return json.dumps({"data": None})
        with open(path, encoding='utf-8') as f:
            return f.read()

    # get arbitrary data
    @app.route("/get_binary_data", methods=["POST"])
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
    @app.route("/tracks/<path:path>")
    def send_track(path):
        file_name = safe_join(project.trackfolder, path)
        range_header = request.headers.get("Range", None)
        if not range_header:
            return send_file(file_name)
        return get_range(file_name, range_header)

    @app.route("/save_state", methods=["POST"])
    def save_data():
        success = True
        try:
            state = request.json
            project.save_state(state)
        except Exception:
            success = False
        return jsonify({"success": success})
    
    return app


def serve_project(project: MDVProject | str , port: int = 5050, open_browser: bool = True ):
    """Serve an MDV project using a lightweight Flask server. 
    Args:
        project (MDVProject | str): The MDV project instance or path to the project directory.
        port (int): The port to run the server on. Defaults to 5050.
        open_browser (bool): Whether to open the browser automatically. Defaults to True.
    """
    if isinstance(project, str):
        project = MDVProject(project)
    app = create_app(project)
    if open_browser:
        import webbrowser
        webbrowser.open(f"http://localhost:{port}/")
    app.run(port=port)


if __name__ == "__main__":
    try:
        if len(sys.argv) > 1:
            path = sys.argv[1]
        else:
            raise ValueError("No project path provided as a command-line argument.")
        if not os.path.exists(path):
            raise FileNotFoundError(f"{path} not found")
        ds_path = os.path.join(path, "datasources.json")
        if not os.path.exists(ds_path):
            raise FileNotFoundError(f"{path} does not contain a valid MDV project.")
        serve_project(MDVProject(path))

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
