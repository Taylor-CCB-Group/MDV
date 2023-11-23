import os
from mdvtools.mdvproject import MDVProject
from .server import add_safe_headers
from flask import Flask,render_template,request,make_response,send_file,Response,jsonify

"""
Make a Flask app, and open all folders in ~/mdv/ as projects that can be served by it.
"""

# todo: make this configurable
project_dir = os.path.join(os.path.expanduser('~'), 'mdv')
# create the directory if it doesn't exist
if not os.path.exists(project_dir):
    os.makedirs(project_dir)

projects = [MDVProject(os.path.join(project_dir, d)) for d in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, d))]

app = Flask(__name__)
app.after_request(add_safe_headers)

for p in projects:
    p.serve(open_browser=False, app=app)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/projects')
def get_projects():
    return jsonify([p.name for p in projects])

if __name__ == '__main__':
    app.run(debug=True, port=5051)