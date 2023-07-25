from flask import Flask,render_template,request,make_response,send_file,Response,jsonify
from flask_socketio import SocketIO
import webbrowser
import mimetypes
import json
import sys
import re
from werkzeug.utils import safe_join
from .websocket import mdv_socketio


def add_safe_headers(resp):
    resp.headers["Cross-Origin-Opener-Policy"]= "same-origin"
    resp.headers["Cross-Origin-Embedder-Policy"]="require-corp"
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

def create_app(project, open_browser=True, port = 5050, websocket=True):
    app=Flask(__name__)
    #add headers to allow web workers
    app.after_request(add_safe_headers)
    if websocket:
      mdv_socketio(app)

    @app.route("/")
    def index():
        return render_template("page.html")
    
    #@app.route("/static/js/<path:path>")
    #def get_js_files(path):
    #    return send_file(safe_join(js_files,path))
    

    @app.route("/<file>.gz")
    def get_binary_file(file):
        # should this now b '.gz'?
        file_name = safe_join(project.dir, file+".gz")
        range_header = request.headers.get('Range', None) 
        return get_range(file_name,range_header)


    @app.route("/<file>.json/")
    def get_json_file(file):
        return send_file(safe_join(project.dir,file+".json"))
      
    #gets the raw byte data and packages it in the correct response
    @app.route('/get_data',methods =["POST"])
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
    @app.route("/images/<path:path>")
    def images(path):
        try:
            return send_file(project.get_image(path))
        except:
            return send_file(safe_join(project.imagefolder, path))

    #All the project's metadata
    @app.route('/get_configs',methods =["GET","POST"])
    def get_configs():
        return json.dumps(project.get_configs())

    #gets a particular view
    @app.route("/get_view",methods = ["POST"])
    def get_view():
        data=request.json
        return json.dumps(project.get_view(data["view"]))
    
    #get any custom row data
    @app.route("/get_row_data",methods=["POST"])
    def get_row_data():
        req=request.json
        try:
            with open( safe_join(project.dir,"rowdata",req["datasource"],f"{req['index']}.json")) as f:
                data = f.read()
        except Exception as e:
            data=json.dumps({"data":None})
        return data
    
    # only the specified region of track files (bam,bigbed,tabix)
    # needs to be returned 
    @app.route("/tracks/<path:path>")
    def send_track(path):
        file_name =safe_join(project.trackfolder,path)
        range_header = request.headers.get('Range', None)  
        if not range_header:
            return send_file(file_name)
        return get_range(file_name,range_header)
    
    @app.route("/save_state",methods = ["POST"])
    def save_data():
        success=True
        try:
            state=request.json
            project.save_state(state)
        except Exception as e:
            success=False

        return jsonify({"success":success})

    if open_browser:
        webbrowser.open(f"http://localhost:{port}")
          
    app.run(port=port)


  

