# project_service.py

from mdvtools.dbutils.dbmodels import db, Project, File
from datetime import datetime
from mdvtools.mdvproject import MDVProject

class ProjectService:
    # list of tuples containing failed project IDs and associated error messages/exceptions
    failed_projects: list[tuple[int, str | Exception]] = []
    @staticmethod
    #  Routes -> /projects
    def get_active_projects():
        try:
            failed_project_ids = {f[0] for f in ProjectService.failed_projects}  
            projects = Project.query.filter(Project.is_deleted == False).all()

            for p in projects:
                print(next(iter(list(MDVProject(p.path).views.values())), None))
            
            # Convert to JSON-ready list of dictionaries
            project_list = [
                {
                    "id": p.id,
                    "name": p.name,
                    "lastModified": p.update_timestamp.strftime('%Y-%m-%d %H:%M:%S'),
                    "thumbnail": next(
                        (v["viewImage"] for v in MDVProject(p.path).views.values() if "viewImage" in v), None
                        ),
                }
                for p in projects if p.id not in failed_project_ids
            ]

            return project_list  #Already formatted for JSON storage
        except Exception as e:
            print(f"Error in dbservice: Error querying active projects: {e}")
            raise

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
            print(f"Error in dbservice: Error getting next project ID: {e}")
            raise

    @staticmethod
    #  Routes -> /create_project
    def add_new_project(path, name='unnamed_project'):
        try:
            new_project = Project(name=name, path=path)
            db.session.add(new_project)
            db.session.commit()
            return new_project
        except Exception as e:
            print(f"Error in dbservice: Error creating project: {e}")
            db.session.rollback()
            raise
        
    @staticmethod
    #  Routes -> /delete_project
    def get_project_by_id(id):
        try:
            return Project.query.get(id)
        except Exception as e:
            print(f"Error in dbservice: Error querying project by id-No project found with id: {e}")
            raise

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
            else:
                print(f"Attempted to soft delete non-existing project with id: {id}")
                return False
        except Exception as e:
            print(f"Error in dbservice: Error soft deleting project: {e}")
            db.session.rollback()
            raise

    @staticmethod
    #  Routes -> /rename_project
    def update_project_name(project_id, new_name):
        try:
            project = Project.query.get(project_id)
            if project:
                project.name = new_name
                # No need to update `update_timestamp` manually
                project.update_timestamp = datetime.now()
                project.accessed_timestamp = datetime.now()
                db.session.commit()
                return True
            return False
        except Exception as e:
            print(f"Error in dbservice: Error renaming project: {e}")
            db.session.rollback()
            raise
    
    @staticmethod
    #  Routes -> /access
    def change_project_access(project_id, new_access_level):
        """Change the access level of a project."""
        try:
            # Retrieve the project by ID
            project = ProjectService.get_project_by_id(project_id)
            if project is None:
                return None, f"Project with ID {project_id} not found", 404

            # Update the access level
            project.access_level = new_access_level
            project.update_timestamp = datetime.now()
            project.accessed_timestamp = datetime.now()
            db.session.commit()
            return project.access_level, "Success", 200

        except Exception as e:
            print(f"Error in dbservice: Unexpected error in change access level: {e}")
            db.session.rollback()
            raise

    @staticmethod
    # Method to set the project's update timestamp
    def set_project_update_timestamp(project_id: str):
        try:
            # Use get_project_by_id to retrieve the project
            project = ProjectService.get_project_by_id(project_id)
            
            if project is None:
                print(f"No project found with ID {project_id}")
                return  # Exit if the project doesn't exist

            # Update the update_timestamp to current datetime
            project.update_timestamp = datetime.now()
            db.session.commit()  # Commit the changes to the database
            print(f"Set project update timestamp for project ID {project_id}")

        except Exception as e:
            print(f"Error in dbservice: Error setting project update timestamp for ID {project_id}: {e}")
            db.session.rollback()
            raise
    
    @staticmethod
    # Method to set the project's accessed timestamp
    def set_project_accessed_timestamp(project_id: str):
        try:
            # Use get_project_by_id to retrieve the project
            project = ProjectService.get_project_by_id(project_id)
            
            if project is None:
                print(f"No project found with ID {project_id}")
                return  # Exit if the project doesn't exist

            # Set the accessed_timestamp to the current datetime
            project.accessed_timestamp = datetime.now()
            db.session.commit()  # Commit the changes to the database
            print(f"Set project accessed timestamp for project ID {project_id}")

        except Exception as e:
            print(f"Error in dbservice: Error setting accessed timestamp for project ID {project_id}: {e}")
            db.session.rollback()
            raise
        
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
                print(f"Added new file to DB: {new_file}")

            # Commit the transaction after adding/updating
            db.session.commit()

        except Exception as e:
            print(f"Error in FileService.add_or_update_file_in_project: Failed to add or update file '{file_name}' for project ID '{project_id}': {str(e)}")
            db.session.rollback()  # Rollback session on error
            raise
        
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

