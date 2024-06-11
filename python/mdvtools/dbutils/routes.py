from flask import request, jsonify, render_template
from mdvtools.dbutils.app import app
from mdvtools.dbutils.dbmodels import db, Project, File
from datetime import datetime
import os

def register_routes(projects):
    
    @app.route('/projects')
    def get_projects():
        return jsonify([p.name for p in projects])
    

def register_global_routes():
    
    @app.route('/')
    def index():
        return render_template('index.html')
    
    @app.route('/upload', methods=['POST'])
    def upload():
        try:
            project_name = request.form.get('project_name')
            if not project_name:
                return jsonify({'error': 'Project name is missing.'}), 400
            
            file = request.files['file']
            if not file:
                return jsonify({'error': 'No file selected.'}), 400

            file_path = os.path.join(app.config['projects_base_dir'], project_name, file.filename)  # type: ignore

            # Check if project exists
            project = Project.query.filter_by(name=project_name).first()
            if not project:
                project = Project(name=project_name)  # type: ignore
                db.session.add(project)
                db.session.commit()  # Commit to get the project.id

            # Check if file with the same name already exists in the project
            existing_file = File.query.filter_by(name=file.filename, project_id=project.id).first()
            if existing_file:
                # Replace the existing file in the file system
                if os.path.exists(existing_file.file_path):
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
                new_file = File(name=file.filename, file_path=file_path, project_id=project.id)  # type: ignore
                
                db.session.add(new_file)
                db.session.commit()

                return jsonify({'message': f'File uploaded successfully under project "{project_name}"'}), 200
        except Exception as e:
            return jsonify({'error in /upload api': str(e)}), 500
