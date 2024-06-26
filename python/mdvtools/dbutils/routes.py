from flask import request, jsonify, render_template
from mdvtools.mdvproject import MDVProject
#from mdvtools.dbutils.app import app
#from mdvtools.dbutils.mdv_server_app import app, db
from mdvtools.dbutils.dbmodels import Project, File
from mdvtools.project_router import ProjectBlueprint
from datetime import datetime
import os


def register_global_routes(app, db, project_dir):
    
    @app.route('/')
    def index():
        return render_template('index.html')
    
    @app.route('/projects')
    def get_projects():
        
        try:
            # Query the database to get all projects that aren't deleted
            print('/projects queried...')
            
            projects = Project.query.filter_by(is_deleted=False).all()
            
            # Return the list of projects with their IDs and names
            return jsonify([{"id": p.id, "name": p.name} for p in projects])
        except Exception as e:
            return jsonify({"status": "error", "message": str(e)}), 500
    
    @app.route("/create_project", methods=["POST"])
    def create_project():
        try:
            print("Creating project")
            
            # Get the next ID from the database
            next_id = db.session.query(db.func.max(Project.id)).scalar()
            if next_id is None:
                next_id = 1
            else:
                next_id += 1
            
            # Create the project directory path
            project_path = os.path.join(project_dir, str(next_id))
            
            # Create and serve the MDVProject
            print("Creating and serving the new project")

            p = MDVProject(project_path)
            p.set_editable(True)
            p.serve(app=app, open_browser=False)
            
            # Create a new Project record in the database with the path
            print("Adding new project to the database")
            new_project = Project.create_project(path=project_path)
            
            return jsonify({"id": new_project.id, "name": new_project.id, "status": "success"})
        except Exception as e:
            return jsonify({"status": "error", "message": str(e)}), 500
        
    
    @app.route("/delete_project/<project_id>", methods=["DELETE"])
    def delete_project(project_id: int):
        """Soft delete a project by setting the deleted flag."""
        try:
            print(f"Deleting project '{project_id}'")
            
            # Check if the project exists in the ProjectBlueprint.blueprints dictionary
            if project_id not in ProjectBlueprint.blueprints:
                return jsonify({"status": "error", "message": f"Project with ID {project_id} not found in ProjectBlueprint.blueprints"}), 404
            else:
                # Find the project by ID and mark it as deleted
                project = Project.query.filter_by(id=project_id).first()
                if project is None:
                    return jsonify({"status": "error", "message": f"Project with ID {project_id} not found in database"}), 404
                
                # Remove the project from the ProjectBlueprint.blueprints dictionary
                del ProjectBlueprint.blueprints[project_id]
                print(f"Removed project '{project_id}' from ProjectBlueprint.blueprints")
                
                
                # Soft delete the project
                project.soft_delete()
                
                return jsonify({"status": "success"})
        except Exception as e:
            return jsonify({"status": "error", "message": str(e)}), 500


    
    
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
