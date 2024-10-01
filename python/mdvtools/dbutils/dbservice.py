# project_service.py

from mdvtools.dbutils.dbmodels import db, Project, File
from datetime import datetime

class ProjectService:
    @staticmethod
    #  Routes -> /projects
    def get_active_projects():
        try:
            return Project.query.filter_by(is_deleted=False).all()
        except Exception as e:
            print(f"Error querying active projects: {e}")
            return []

    @staticmethod
    #  Routes -> /create_project
    def get_next_project_id():
        try:
            next_id = db.session.query(db.func.max(Project.id)).scalar()
            if next_id is None:
                next_id = 1
            else:
                next_id += 1
            return next_id
        except Exception as e:
            print(f"Error getting next project ID: {e}")
            return None

    @staticmethod
    #  Routes -> /create_project
    def add_new_project(path, name='unnamed_project'):
        try:
            new_project = Project(name=name, path=path)
            db.session.add(new_project)
            db.session.commit()
            return new_project
        except Exception as e:
            print(f"Error creating project: {e}")
            return None
        
    @staticmethod
    #  Routes -> /delete_project
    def get_project_by_id(id):
        try:
            return Project.query.get(id)
        except Exception as e:
            print(f"Error querying project by id: {e}")
            return None

    @staticmethod
    #  Routes -> /delete_project
    def soft_delete_project(id):
        try:
            project = Project.query.get(id)
            if project:
                project.is_deleted = True
                project.deleted_timestamp = datetime.now()
                db.session.commit()
                return True
            return False
        except Exception as e:
            print(f"Error soft deleting project: {e}")
            return False
        
class FileService:
    @staticmethod
    def add_or_update_file_in_project(file_name, file_path, project_id):
        """Adds a new file or updates an existing file in the database."""
        try:
            existing_file = FileService.get_file_by_path_and_project(file_path, project_id)

            if existing_file:
                # Update the file details if they are different
                if existing_file.name != file_name:
                    existing_file.name = file_name
                    existing_file.update_timestamp = datetime.now()
                    print(f"Updated file name in DB: {existing_file}")
                db.session.commit()
                return existing_file
            else:
                # Add new file to the database
                new_file = File(
                    name=file_name,
                    file_path=file_path,
                    project_id=project_id,
                    upload_timestamp=datetime.now(),
                    update_timestamp=datetime.now()
                )
                db.session.add(new_file)
                db.session.commit()
                print(f"Added new file to DB: {new_file}")
                return new_file
        except Exception as e:
            print(f"Error adding or updating file in project: {e}")
            db.session.rollback()  # Rollback session on error
            return None
        
    @staticmethod
    def get_file_by_path_and_project(file_path, project_id):
        """Fetch a file by its path and project ID."""
        try:
            return File.query.filter_by(file_path=file_path, project_id=project_id).first()
        except Exception as e:
            print(f"Error retrieving file by path '{file_path}' and project ID {project_id}: {e}")
            return None
    
    @staticmethod
    def file_exists_in_project(file_path, project_id):
        """Utility function to check if a file exists in the files table."""
        try:
            return File.query.filter_by(
                file_path=file_path, project_id=project_id
            ).first() is not None
        except Exception as e:
            print(f"Error checking file existence: {e}")
            return False

    @staticmethod
    def get_files_by_project(project_id):
        try:
            return File.query.filter_by(project_id=project_id).all()
        except Exception as e:
            print(f"Error querying files for project ID {project_id}: {e}")
            return []

    @staticmethod
    def delete_files_by_project(project_id):
        try:
            files = File.query.filter_by(project_id=project_id).all()
            for file in files:
                db.session.delete(file)
            db.session.commit()
            return True
        except Exception as e:
            print(f"Error deleting files for project ID {project_id}: {e}")
            db.session.rollback()  # Rollback session on error
            return False

    @staticmethod
    def update_file_timestamp(file_id):
        try:
            file = File.query.get(file_id)
            if file:
                file.update_timestamp = datetime.now()
                db.session.commit()
                return True
            return False
        except Exception as e:
            print(f"Error updating timestamp for file ID {file_id}: {e}")
            db.session.rollback()  # Rollback session on error
            return False

