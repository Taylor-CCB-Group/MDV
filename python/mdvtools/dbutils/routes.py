from flask import request, jsonify, render_template
from mdvtools.dbutils.app import app
from mdvtools.dbutils.dbmodels import db, Project, File
from datetime import datetime
import os

def register_routes(projects):
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