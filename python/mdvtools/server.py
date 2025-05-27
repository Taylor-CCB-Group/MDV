from flask import (
    Flask,
    # Blueprint,
    render_template,
    request,
    make_response,
    send_file,
    Response,
    jsonify,
    current_app
)
import webbrowser
import mimetypes
import json
import sys
import re
from mdvtools.llm.code_manipulation import parse_view_name
from werkzeug.security import safe_join
from mdvtools.websocket import mdv_socketio
from mdvtools.mdvproject import MDVProject
from mdvtools.project_router import (
    ProjectBlueprint as Blueprint,
    SingleProjectShim
)
import os
import pandas as pd
from datetime import datetime
from typing import Optional
import threading
import scanpy as sc
from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.llm.chat_protocol import ProjectChatProtocol, ProjectChat
from mdvtools.llm.chat_server_extension import chat_extension
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
    # todo: another way of specifying that we are using chat... also, use socket for other things.
    websocket=True, # todo - pass something back to client in `state.json` that indicates whether this is enabled.
    use_reloader=False,
    app: Optional[Flask] = None,
    # would like to redesign this a bit... this relates mostly to the ProjectBlueprint "_v2" thing
    # and also some BS about adding headers to cope with server misconfiguration (as I recall).
    backend_db=False,
):
    if app is None:
        route = ""
        # route = "/project/" + project.name # for testing new API with simple app...
        app = Flask(__name__)
        print(f"created Flask {app}")
        # add headers to allow web workers
        app.after_request(add_safe_headers)
        project_bp = SingleProjectShim(app)
        multi_project = False
        if websocket:
            # reviewing this... thinking about hooking up to ProjectChat logger...
            #! nb - we're in 'single project' mode here.
            # thinking about syncronising list of views via SocketIO rather than polling... 
            # maybe via a smaller PR where I better figure out clean socket implementation.
            # maybe we have some abstraction around how we create the Flask instance.
            mdv_socketio(app)
    else:
        ## nb - previous use of flask.Blueprint was not allowing new projects at runtime
        ## we substitute this with our own ProjectBlueprint class, which is a drop-in replacement
        ## but we should add more tests to ensure it behaves as expected...
        # add routes for this project to existing app
        # set the route prefix to the project name, derived from the dir name.
        # this is to allow multiple projects to be served from the same server.
        multi_project = True
        route = "/project/" + project.id + "/"
        

        if backend_db:
            from mdvtools.project_router import ProjectBlueprint_v2 as Blueprint_v2
            print("backend_db is True")
            project_bp = Blueprint_v2(project.id, __name__, url_prefix=route)
        else:
            project_bp = Blueprint(project.id, __name__, url_prefix=route)

    # if route in routes:
    #     raise Exception(
    #         "Route already exists - can't have two projects with the same name"
    #     )
    routes.add(route)

    if websocket:
        """
        current prototype using 'websocket' as a flag that we interpret as 'enable chatMDV' for now.
        we might avoid having this logic in here, or at least change that semantic.
        rather than pass in `app: Optional[Flask]`, we could pass in a configuration object
        that includes that app and any extensions we want to use.
        for now, the code has been rearranged, but the logic for calling it is similar...
        """
        chat_extension.register_routes(project, project_bp)


    @project_bp.route("/")
    def project_index():
        print("recieved request to project_index")
        # the backend page currently needs to be different to workaround a server config issue
        # some requests were being downgraded to http, which caused problems with the backend
        # but if we always add the header it messes up localhost development.
        # todo if necessary, apply equivalent change to index.html / any other pages we might have
        return render_template("page.html", route=route, backend=backend_db)

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
        # print(f"get_json_file: '{path}' for project {project.id}")

        
        if path is None or not os.path.exists(path):
            return "File not found", 404
        if file == "state":
            with open(path) as f:
                try:
                    state = json.load(f)
                    # do we want this to always be true/not a flag we pass?
                    state["websocket"] = websocket
                    # in future, we could iterate over a list of extensions.
                    chat_extension.mutate_state_json(state, project)
                    return state
                except Exception as e:
                    return f"Problem parsing state file: {e}", 500
        return _send_file(path)

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

    @project_bp.route("/save_state", access_level='editable', methods=["POST"])
    def save_data():
        success = True
        try:
            state = request.json
            project.save_state(state)
        except Exception:
            success = False

        return jsonify({"success": success})

    # Utility Functions
    def create_temp_folder(base_path):
        """
        Create a temporary folder with a timestamp-based name.

        Creates a new directory within the specified base path using a timestamp-based naming
        convention for temporary AnnData file storage. The folder name follows the pattern
        'temp_anndata_YYYYMMDD_HHMMSS'.

        Args:
            base_path (str): The parent directory where the temporary folder will be created

        Returns:
            str: Absolute path to the created temporary folder

        Notes:
            - The folder is created with exist_ok=True to handle potential race conditions
            - Timestamp format used: YYYYMMDD_HHMMSS (e.g., temp_anndata_20250129_143022)
            - This function is typically used in conjunction with cleanup_folder() for
            temporary file handling during AnnData uploads

        Example:
            >>> temp_path = create_temp_folder('/path/to/project')
            >>> print(temp_path)
            '/path/to/project/temp_anndata_20250129_143022'
        """
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        temp_folder_name = f"temp_anndata_{timestamp}"
        temp_folder_path = os.path.join(base_path, temp_folder_name)
        os.makedirs(temp_folder_path, exist_ok=True)
        return temp_folder_path

    def cleanup_folder(folder_path):
        """
        Recursively delete a folder and all its contents.

        Safely removes all files within the specified folder and then removes the folder
        itself. Any errors during the cleanup process are logged but do not raise exceptions
        to the caller.

        Args:
            folder_path (str): Absolute path to the folder that needs to be cleaned up

        Notes:
            - First removes all files within the folder, then removes the empty folder
            - Handles non-existent paths safely
            - Logs errors to current_app.logger but does not propagate exceptions
            - Only removes files (not subdirectories) within the specified folder
            - Typically used to clean up temporary folders created by create_temp_folder()

        Example:
            >>> cleanup_folder('/path/to/temp/folder')
            # Folder and its contents are deleted if they exist
            # Any errors are logged but not raised
        """
        try:
            if os.path.exists(folder_path):
                for file in os.listdir(folder_path):
                    file_path = os.path.join(folder_path, file)
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                os.rmdir(folder_path)
        except Exception as e:
            current_app.logger.error(f"Error cleaning up folder {folder_path}: {e}")

    @project_bp.route("/add_anndata", access_level='editable', methods=["POST"])
    def add_anndata():
        """
        Upload and process an AnnData file (.h5ad) for the project.

        This endpoint handles the upload ofa single file in AnnData format,
        converts it to MDV format, and saves it to the project directory. If a file already exists,
        it will be saved to a temporary location that expires after 5 minutes.

        HTTP Methods:
            POST

        Request Parameters:
            file (FileStorage): The .h5ad file to be uploaded (multipart/form-data)

        Returns:
            tuple: A tuple containing:
                - JSON response with status and message
                - HTTP status code

            Success Response (200):
                {
                    'status': 'success',
                    'message': 'Anndata added successfully'
                }

            Error Responses:
                400:
                    - No file in request:
                        {
                            'status': 'error',
                            'message': 'No file part in the request'
                        }
                    - No file selected:
                        {
                            'status': 'error',
                            'message': 'No file selected'
                        }
                    - Invalid file type:
                        {
                            'status': 'error',
                            'message': 'Invalid file type. Only .h5ad files are allowed'
                        }
                409:
                    - File exists conflict:
                        {
                            'status': 'conflict',
                            'message': 'File already exists',
                            'temp_folder': '/path/to/temp/folder'
                        }

        Notes:
            - The function creates a temporary folder for conflict resolution that is automatically
            cleaned up after 5 minutes (300 seconds)
            - Only .h5ad file extensions are accepted
            - Successful upload will trigger conversion from Scanpy AnnData to MDV format
        """
        try:
            if 'file' not in request.files:
                return jsonify({'status': 'error', 'message': 'No file part in the request'}), 400

            file = request.files['file']
            if not file or not file.filename:
                return jsonify({'status': 'error', 'message': 'No file selected'}), 400

            if not file.filename.endswith('.h5ad'):
                return jsonify({'status': 'error', 'message': 'Invalid file type. Only .h5ad files are allowed'}), 400

            target_path = os.path.join(project.dir, "anndata.h5ad")
            if os.path.exists(target_path):
                temp_folder = create_temp_folder(project.dir)
                temp_path = os.path.join(temp_folder, "anndata.h5ad")
                file.save(temp_path)

                threading.Timer(300, cleanup_folder, args=[temp_folder]).start()
                return jsonify({
                    'status': 'conflict',
                    'message': 'File already exists',
                    'temp_folder': temp_folder
                }), 409

            # Save file and process
            file.save(target_path)
            anndata = sc.read(target_path)
            convert_scanpy_to_mdv(project.dir, anndata)
            return jsonify({'status': 'success', 'message': 'Anndata added successfully'}), 200

        except Exception as e:
            current_app.logger.error(f"Unexpected error: {str(e)}")
            return jsonify({'status': 'error', 'message': str(e)}), 500

    @project_bp.route("/combine_anndata", access_level='editable', methods=["PATCH"])
    def combine_anndata():
        """
        Combine a temporary AnnData file with an existing project file.

        This endpoint handles the resolution of file conflicts when uploading AnnData files.
        It can either combine the temporary file with the existing project file (adding a label prefix
        to distinguish the data) or cancel the operation and clean up the temporary files.

        HTTP Methods:
            PATCH

        Request Parameters:
            temp_folder (str): Path to the temporary folder containing the AnnData file
            combine (bool): Whether to combine the files (true) or cancel the operation (false)
            label (str): Prefix label to distinguish the new data when combining files

        Returns:
            tuple: A tuple containing:
                - JSON response with status and message
                - HTTP status code

            Success Responses:
                200:
                    - Successful merge:
                        {
                            'status': 'success',
                            'message': 'File merged successfully'
                        }
                    - Operation cancelled:
                        {
                            'status': 'success',
                            'message': 'Operation cancelled'
                        }

            Error Responses:
                400:
                    - Missing label:
                        {
                            'status': 'error',
                            'message': 'Label field not found'
                        }
                    - Missing temporary file:
                        {
                            'status': 'error',
                            'message': 'Temporary file not found'
                        }
                408:
                    - Timeout:
                        {
                            'status': 'error',
                            'message': 'Request timed out, please try uploading again'
                        }
                500:
                    - Label generation error:
                        {
                            'status': 'error',
                            'message': '<specific ValueError message>'
                        }
                    - Other errors:
                        {
                            'status': 'error',
                            'message': '<specific error message>'
                        }

        Notes:
            - This endpoint is typically called after a conflict is detected in the /add_anndata endpoint
            - The temporary folder is cleaned up regardless of whether the operation succeeds or fails
            - When combining files, the new data is prefixed with the provided label followed by an underscore
            - The combined file maintains the original filename 'anndata.h5ad' in the project directory
        """
        try:
            temp_folder = request.form.get('temp_folder')
            combine = request.form.get('combine') == 'true'
            label = request.form.get('label')

            if not temp_folder or not os.path.exists(temp_folder):
                return jsonify({'status': 'error', 'message': 'Request timed out, please try uploading again'}), 408

            if not combine:
                cleanup_folder(temp_folder)
                return jsonify({'status': 'success', 'message': 'Operation cancelled'}), 200
            
            if not label:
                cleanup_folder(temp_folder)
                return jsonify({'status': 'error', 'message': 'Label field not found'}), 400

            temp_path = os.path.join(temp_folder, "anndata.h5ad")
            if not os.path.exists(temp_path):
                return jsonify({'status': 'error', 'message': 'Temporary file not found'}), 400


            # Process and combine the file
            new_anndata = sc.read(temp_path)
            convert_scanpy_to_mdv(project.dir, new_anndata, delete_existing=False, label=f"{label}_")
            new_anndata.write(os.path.join(project.dir, "anndata.h5ad")) # type: ignore - str output from os.path.join should be reliably PathLike
            cleanup_folder(temp_folder)

            return jsonify({'status': 'success', 'message': 'File merged successfully'}), 200

        except ValueError as ve:
            current_app.logger.error(f"Label generation error: {str(ve)}")
            return jsonify({'status': 'error', 'message': str(ve)}), 500
        except Exception as e:
            current_app.logger.error(f"Unexpected error: {str(e)}")
            return jsonify({'status': 'error', 'message': str(e)}), 500
    
    @project_bp.route("/add_or_update_image_datasource", access_level='editable', methods=["POST"])
    def add_or_update_image_datasource():
        try:
            # Check if request has a file part
            if 'file' not in request.files:
                return "No file part in the request", 400

            # Get the file from the request
            file = request.files['file']
            
            # Get the text fields from the request form
            datasource_name = request.form.get('datasourceName') # ""
            tiff_metadata = request.form.get('tiffMetadata')
            # Validate the presence of required fields
            if not file or not tiff_metadata:
                return "Missing file or tiffMetadata", 400
            # If tiff_metadata is sent as JSON string, deserialize it
            try:
                tiff_metadata = json.loads(tiff_metadata)
            except Exception as e:
                return jsonify({"status": "error", "message": f"Invalid JSON format for tiffMetadata: {e}"}), 400
            
            # Call the method to add or update the image datasource
            view_name = project.add_or_update_image_datasource(tiff_metadata, datasource_name, file)
            
            # If no exception is raised, the operation was successful. let the client know which view will show the image.
            print(f">>> notify client that image datasource updated and file uploaded successfully, view: {view_name}")
            return jsonify({"status": "success", "message": "Image datasource updated and file uploaded successfully", "view": view_name}), 200

        except Exception as e:
            return jsonify({"status": "error", "message": str(e)}), 500


    @project_bp.route("/add_datasource", access_level='editable', methods=["POST"])
    def add_datasource():
        # we shouldn't be passing "backend" in request.form, the logic should only be on server
        #if backend:
        #    response = add_datasource_backend(project)
        #    return response

        if (
            "permission" not in project.state
            or not project.state["permission"] == "edit"
        ):
            return "Project is read-only", 400
        success = True
        try:
            name = request.form["name"]

            print("In server.py add_datasource")

            if not name:
                return "Request must contain 'name'", 400
            # xxx - not how column metadata should be passed, todo fix
            # cols = (
            #     request.form["columns"].split(",")
            #     if "columns" in request.form
            #     else None
            # )

            # I'm not sure we really want to add to default view by default - could mess up existing views in a project with multiple datasources
            # but probably ok for now (famous last words)
            view = request.form["view"] if "view" in request.form else "default"
            # replace = True if "replace" in request.form else False
            replace = False
            if not replace and name in [ds["name"] for ds in project.datasources]:
                return (
                    f"Datasource '{name}' already exists, and 'replace' was not set in request",
                    400,
                )
            if "file" not in request.files:
                return "No 'file' provided in request form data", 400
            file = request.files["file"]
            supplied_only = True if "supplied_only" in request.form else False
            if not file or file.mimetype != "text/csv":
                return "File must be a CSV", 400
            file.seek(0)
            # will this work? can we return progress to the client?
            file_name = project.dir + "/table.csv"
            file.save(file_name)
            #df = pd.read_csv(file.stream)
            df = pd.read_csv(file_name)
            print("In server.py add_datasource- df created")
            print("df is ready, calling project.add_datasource")
            project.add_datasource(
                #project.id,
                name,
                df,
                # cols,
                add_to_view=view,
                supplied_columns_only=supplied_only,
                replace_data=replace
                )
            print("added df - project.add_datasource completed")
        except Exception as e:
            # success = False
            return str(e), 400

        metadata = project.get_datasource_metadata(name)
        return jsonify({"success": success, "metadata": metadata})

    if open_browser:
        webbrowser.open(f"http://localhost:{port}/{route}")

    if multi_project:
        assert(isinstance(app, Flask))
        if route in app.blueprints:
            print(f"there is already a blueprint at {route}")
        print(f"Adding project {project.id} to existing app")
        ## nb - uncomment this if not using ProjectBlueprint refactor...
        # app.register_blueprint(project_bp)
    else:
        # user_reloader=False, allows the server to work within jupyter

        # app.run(host="0.0.0.0", port=port, debug=True, use_reloader=use_reloader)
        ## todo - gevent for mdvlite / non-optional dependency
        from gevent.pywsgi import WSGIServer
        http_server = WSGIServer(("127.0.0.1", port), app)
        http_server.serve_forever()


