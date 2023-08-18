from flask import Flask,render_template,request,make_response,send_file,Response,jsonify
from flask_socketio import SocketIO
import webbrowser
import mimetypes
import json
import sys
import re
from werkzeug.utils import safe_join
from .websocket import mdv_socketio


routes = set()

def add_safe_headers(resp):
    resp.headers["Cross-Origin-Opener-Policy"]= "same-origin"
    resp.headers["Cross-Origin-Embedder-Policy"]="require-corp"
    resp.headers["Access-Control-Allow-Origin"]="*"
    return resp



def get_range(file_name,range_header):
    file =open(file_name,"rb")
    size = sys.getsizeof(file_name)
    byte1, byte2 = 0, None

    m = re.search('(\d+)-(\d*)', range_header)
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
    rv = Response(data,
                206,
                mimetype=mimetypes.guess_type(file_name)[0],
                direct_passthrough=True)
    rv.headers.add('Content-Range', 'bytes {0}-{1}/{2}'.format(byte1, byte1 + length - 1, size))
    rv.headers.add('Accept-Ranges', 'bytes')
    file.close()
    return rv

def create_app(project, open_browser = True, port = 5050, websocket = False, app:Flask = None):
    if app is None:
        route = ""
        # route = "/project/" + project.name # for testing new API with simple app...
        app=Flask(__name__)
        #add headers to allow web workers
        app.after_request(add_safe_headers)
        # nb, may make this default to False again soon.
        ### 'MEW' in Unity is using IWeb PostMessage, not WebSockets.
        ### but this will be used for local testing in short-term, and potentially other things later.
        if websocket:
            mdv_socketio(app)
    else:
        # add routes for this project to existing app
        # set the route prefix to the project name, derived from the dir name.
        # this is to allow multiple projects to be served from the same server.
        route = "/project/" + project.name
    if route in routes:
        raise Exception("Route already exists - can't have two projects with the same name")
    routes.add(route)

    @app.route(f"{route}/")
    def index():
        # `_mdvInit('{{route}}')` in template...
        return render_template("page.html", route=route)
    
    #@app.route(f"{route}/static/js/<path:path>")
    #def get_js_files(path):
    #    return send_file(safe_join(js_files,path))

    @app.route(f"{route}/static/<path:path>")
    def get_static_files(path):
        return send_file(safe_join('./static/', path))
    

    @app.route(f"{route}/<file>.b")
    def get_binary_file(file):
        # should this now b '.gz'?
        file_name = safe_join(project.dir, file+".b")
        range_header = request.headers.get('Range', None) 
        return get_range(file_name,range_header)
    
    # duplicate of above, but for .gz files in case that's needed.
    # (there was some reason for changing to this, but I can't fully remember the status
    # so maybe better to support both for now)
    @app.route(f"{route}/<file>.gz")
    def get_binary_file_gz(file):
        file_name = safe_join(project.dir, file+".gz")
        range_header = request.headers.get('Range', None) 
        return get_range(file_name,range_header)


    @app.route(f"{route}/<file>.json")
    def get_json_file(file):
        return send_file(safe_join(project.dir,file+".json"))
      
    #gets the raw byte data and packages it in the correct response
    @app.route(f"{route}/get_data", methods=["POST"])
    def get_data():
        try:
            data = request.json
            bytes_ = project.get_byte_data(data["columns"],data["data_source"])
            response=make_response(bytes_)
            response.headers.set('Content-Type', 'application/octet-stream')
            return response
        except Exception as e:
            print(e) 
            return "Problem handling request",400
    
    #images contained in the project
    @app.route(f"{route}/images/<path:path>")
    def images(path):
        try:
            return send_file(project.get_image(path))
        except:
            return send_file(safe_join(project.imagefolder, path))

    #All the project's metadata
    @app.route(f"{route}/get_configs", methods=["GET","POST"])
    def get_configs():
        return json.dumps(project.get_configs())

    #gets a particular view
    @app.route(f"{route}/get_view", methods=["POST"])
    def get_view():
        data=request.json
        return json.dumps(project.get_view(data["view"]))
    
    #get any custom row data
    @app.route(f"{route}/get_row_data", methods=["POST"])
    def get_row_data():
        req=request.json
        try:
            with open( safe_join(project.dir,"rowdata",req["datasource"],f"{req['index']}.json")) as f:
                data = f.read()
        except Exception as e:
            data=json.dumps({"data":None})
        return data
    
    #get arbitrary data
    @app.route(f"{route}/get_binary_data", methods=["POST"])
    def get_binary_data():
        req=request.json
        try:
            with open( safe_join(project.dir,"binarydata",req["datasource"],f"{req['name']}.b"),"rb") as f:
                data = f.read()
        except Exception as e:
            data=None
        return data

    # only the specified region of track files (bam,bigbed,tabix)
    # needs to be returned 
    @app.route(f"{route}/tracks/<path:path>")
    def send_track(path):
        file_name =safe_join(project.trackfolder,path)
        range_header = request.headers.get('Range', None)  
        if not range_header:
            return send_file(file_name)
        return get_range(file_name,range_header)
    
    @app.route(f"{route}/save_state", methods = ["POST"])
    def save_data():
        success=True
        try:
            state=request.json
            project.save_state(state)
        except Exception as e:
            success=False

        return jsonify({"success":success})

    if open_browser:
        webbrowser.open(f"http://localhost:{port}/{route}")
          
    app.run(port=port)


  

