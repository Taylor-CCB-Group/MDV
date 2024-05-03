import os
from mdvtools.mdvproject import MDVProject
from ..server import add_safe_headers
from flask import Flask,render_template,request,jsonify
from flask_sqlalchemy import SQLAlchemy
from datetime import datetime
import json


"""
Make a Flask app, and open all folders in ~/mdv/ as projects that can be served by it.
"""

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://postgres@localhost/mydatabase'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['UPLOAD_FOLDER'] = 'uploads'

db = SQLAlchemy(app)

class Project(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(255), nullable=False, unique=True)
    files = db.relationship('File', backref='project', lazy=True)

class File(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(255), nullable=False)
    file_path = db.Column(db.String(255), nullable=False)
    upload_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now)
    update_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    project_id = db.Column(db.Integer, db.ForeignKey('project.id'), nullable=False)

class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(80), unique=True, nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=False)


# todo: make this configurable
project_dir = os.path.join(os.path.expanduser('~'), 'mdv')
# create the directory if it doesn't exist
if not os.path.exists(project_dir):
    os.makedirs(project_dir)

# list of projects to serve needs to come from db
projects = [MDVProject(os.path.join(project_dir, d)) for d in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, d))]


# add other projects as listed in config file (maybe later a sqlite db or something)
# nothing more permanent than the provisional... just want a quick way of including projects on an external volume.
# manually editted for now, and not watched for changes.
config_file = os.path.join(project_dir, 'config.json')
if os.path.exists(config_file):
    import json
    with open(config_file, 'r') as f:
        config = json.load(f)
    for d in config['projects']:
        print(f"adding project '{d}' from config file")
        try:
            projects.append(MDVProject(d))
        except Exception as e:
            print(f"error '{e}' adding project '{d}' from config file")
else:
    with open(config_file, 'w') as f:
        json.dump({'projects': []}, f)


@app.route('/')
def index():
    # todo: figure out what to do here / how to configure routes
    return render_template('index.html')

@app.route('/projects')
def get_projects():
    return jsonify([p.name for p in projects])

@app.route('/upload', methods=['POST'])
def upload():
    try:
        project_name = request.form.get('project_name')
        if not project_name:
            return jsonify({'error': 'Project name is missing.'}), 400
        
        file = request.files['file']
        if not file:
            return jsonify({'error': 'No file selected.'}), 400

        file_path = os.path.join(app.config['UPLOAD_FOLDER'], file.filename) # type: ignore

        # Check if project exists
        project = Project.query.filter_by(name=project_name).first()
        if not project:
            project = Project(name=project_name) # type: ignore
            db.session.add(project)

        # Check if file with same name already exists in the project
        existing_file = File.query.filter_by(name=file.filename, project_id=project.id).first()
        if existing_file:
            # Replace the existing file in the file system
            os.remove(existing_file.file_path)

            # Update the database entry with new file path and update timestamp
            existing_file.file_path = file_path
            existing_file.update_timestamp = datetime.now()

            # Save the file to the uploads directory
            file.save(file_path)

            db.session.commit()

            return jsonify({'message': f'File "{file.filename}" under project "{project_name}" exists already. File has been replaced.'}), 200
        else:
            # Save the file to the uploads directory
            file.save(file_path)

            # Create a new file entry
            new_file = File(name=file.filename, file_path=file_path, project=project) # type: ignore
            db.session.add(new_file)
            db.session.commit()

            return jsonify({'message': f'File uploaded successfully under project "{project_name}"'}), 200
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/delete_file', methods=['DELETE'])
def delete_file():
    try:
        file_id = request.args.get('file_id')
        if not file_id:
            return jsonify({'error': 'File ID is missing.'}), 400

        # Check if the file exists
        file = File.query.get(file_id)
        if not file:
            return jsonify({'error': 'File not found.'}), 404

        # Delete the file from the uploads directory
        if os.path.exists(file.file_path):
            os.remove(file.file_path)
        
        # Delete the file entry from the database
        db.session.delete(file)
        db.session.commit()
        return jsonify({'message': 'File deleted successfully.'}), 200
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    
    app.after_request(add_safe_headers)

    with app.app_context():
        db.create_all()

        # Fetch project names from the database
        #projects= Project.query.with_entities(Project.name).all()


        # Extract project names from the list of tuples
        #projects = [project.name for project in projects_entities]

        #print(projects)

    for p in projects:
        try:
            p.serve(open_browser=False, app=app)
        except Exception as e:
            print(f"error '{e}' serving {p.name}...")

    app.run(debug=True, port=5051)