# project_service.py

from mdvtools.dbutils.dbmodels import db, Project, File, User, UserProject
from datetime import datetime
from mdvtools.mdvproject import MDVProject
from typing import Optional
from mdvtools.logging_config import get_logger

logger = get_logger(__name__)

class ProjectService:
    # list of tuples containing failed project IDs and associated error messages/exceptions
    failed_projects: list[tuple[int, str | Exception]] = []
    @staticmethod
    #  Routes -> /projects
    def get_active_projects():
        try:
            failed_project_ids = {f[0] for f in ProjectService.failed_projects}  
            projects = Project.query.filter(~Project.is_deleted).all()


            # Convert to JSON-ready list of dictionaries
            def get_project_thumbnail(project_path):
                """Extract the first available viewImage from a project's views."""
                try:
                    mdv_project = MDVProject(project_path)
                    return next((v["viewImage"] for v in mdv_project.views.values() if "viewImage" in v), None)
                except Exception as e:
                    logger.exception(f"Error extracting thumbnail for project at {project_path}: {e}")
                    return None
            
            def get_readme_file(project_path):
                try:
                    mdv_project = MDVProject(project_path)
                    readme_file = mdv_project.readme
                    return readme_file
                except Exception as e:
                    logger.exception(f"Error getting readme file for project at {project_path}: {e}")
                    return None

            project_list = [
                {
                    "id": p.id,
                    "name": p.name,
                    "lastModified": p.update_timestamp.strftime('%Y-%m-%d %H:%M:%S'),
                    "thumbnail": get_project_thumbnail(p.path),
                    "readme": get_readme_file(p.path),
                }
                for p in projects if p.id not in failed_project_ids
            ]

            return project_list
        except Exception as e:
            logger.exception(f"Error in dbservice: Error querying active projects: {e}")
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
            logger.exception(f"Error in dbservice: Error getting next project ID: {e}")
            raise

    @staticmethod
    #  Routes -> /create_project
    def add_new_project(path, name='unnamed_project'):
        try:
            # new_project = Project(name=name, path=path)
            new_project = Project()
            new_project.name = name
            new_project.path = path
            db.session.add(new_project)
            db.session.commit()
            return new_project
        except Exception as e:
            logger.exception(f"Error in dbservice: Error creating project: {e}")
            db.session.rollback()
            raise
        
    @staticmethod
    #  Routes -> /delete_project
    def get_project_by_id(id):
        try:
            return Project.query.get(id)
        except Exception as e:
            logger.exception(f"Error in dbservice: Error querying project by id-No project found with id: {e}")
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
                logger.info(f"Attempted to soft delete non-existing project with id: {id}")
                return False
        except Exception as e:
            logger.exception(f"Error in dbservice: Error soft deleting project: {e}")
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
            logger.exception(f"Error in dbservice: Error renaming project: {e}")
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
            logger.exception(f"Error in dbservice: Unexpected error in change access level: {e}")
            db.session.rollback()
            raise

    @staticmethod
    # Method to set the project's update timestamp
    def set_project_update_timestamp(project_id: str):
        try:
            # Use get_project_by_id to retrieve the project
            project = ProjectService.get_project_by_id(project_id)
            
            if project is None:
                logger.info(f"No project found with ID {project_id}")
                return  # Exit if the project doesn't exist

            # Update the update_timestamp to current datetime
            project.update_timestamp = datetime.now()
            db.session.commit()  # Commit the changes to the database
            logger.info(f"Set project update timestamp for project ID {project_id}")

        except Exception as e:
            logger.exception(f"Error in dbservice: Error setting project update timestamp for ID {project_id}: {e}")
            db.session.rollback()
            raise
    
    @staticmethod
    # Method to set the project's accessed timestamp
    def set_project_accessed_timestamp(project_id: str):
        try:
            # Use get_project_by_id to retrieve the project
            project = ProjectService.get_project_by_id(project_id)
            
            if project is None:
                logger.info(f"No project found with ID {project_id}")
                return  # Exit if the project doesn't exist

            # Set the accessed_timestamp to the current datetime
            project.accessed_timestamp = datetime.now()
            db.session.commit()  # Commit the changes to the database
            logger.info(f"Set project accessed timestamp for project ID {project_id}")

        except Exception as e:
            logger.exception(f"Error in dbservice: Error setting accessed timestamp for project ID {project_id}: {e}")
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
                    logger.info(f"Updated file name in DB: {existing_file}")
            else:
                # Add new file to the database
                # new_file = File(
                #     name=file_name,
                #     file_path=file_path,
                #     project_id=project_id,
                #     upload_timestamp=datetime.now(),
                #     update_timestamp=datetime.now()
                # )
                new_file = File()
                new_file.name = file_name
                new_file.file_path = file_path
                new_file.project_id = project_id
                new_file.upload_timestamp = datetime.now()
                new_file.update_timestamp = datetime.now()
                db.session.add(new_file)
                logger.info(f"Added new file to DB: {new_file}")

            # Commit the transaction after adding/updating
            db.session.commit()

        except Exception as e:
            logger.exception(f"Error in FileService.add_or_update_file_in_project: Failed to add or update file '{file_name}' for project ID '{project_id}': {str(e)}")
            db.session.rollback()  # Rollback session on error
            raise
        
    @staticmethod
    def get_file_by_path_and_project(file_path, project_id):
        """Fetch a file by its path and project ID."""
        try:
            return File.query.filter_by(file_path=file_path, project_id=project_id).first()
        except Exception as e:
            logger.exception(f"Error retrieving file by path '{file_path}' and project ID {project_id}: {e}")
            return None
    
    @staticmethod
    def file_exists_in_project(file_path, project_id):
        """Utility function to check if a file exists in the files table."""
        try:
            return File.query.filter_by(
                file_path=file_path, project_id=project_id
            ).first() is not None
        except Exception as e:
            logger.exception(f"Error checking file existence: {e}")
            return False

    @staticmethod
    def get_files_by_project(project_id):
        try:
            return File.query.filter_by(project_id=project_id).all()
        except Exception as e:
            logger.exception(f"Error querying files for project ID {project_id}: {e}")
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
            logger.exception(f"Error deleting files for project ID {project_id}: {e}")
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
            logger.exception(f"Error updating timestamp for file ID {file_id}: {e}")
            db.session.rollback()  # Rollback session on error
            return False


class UserService:
    @staticmethod
    def add_or_update_user(email: str, auth_id: Optional[str] = None, first_name: Optional[str] = None, last_name: Optional[str] = None, institution: Optional[str] = None):
        """
        Adds a new user or updates an existing user based on the provided email.

        :param email: User's email address (mandatory).
        :param auth_id: User's Auth ID (optional).
        :param first_name: User's first name (optional).
        :param last_name: User's last name (optional).
        :param institution: User's institution or association (optional).
        :return: The created or updated User object.
        """
        try:
            if not email:
                raise ValueError("Email is required.")

            user = User.query.filter_by(email=email).first()

            if user:
                # Update existing user with provided non-empty fields
                if auth_id:
                    user.auth_id = auth_id
                if first_name:
                    user.first_name = first_name
                if last_name:
                    user.last_name = last_name
                if institution:
                    user.institution = institution
                db.session.commit()
                return user
            else:
                # Create new user with provided fields
                # new_user = User(
                #     email=email,
                #     auth_id=auth_id or '',
                #     first_name=first_name or '',
                #     last_name=last_name or '',
                #     institution=institution,
                #     confirmed_at=datetime.utcnow(),
                #     is_active=True,
                #     password='',  # Set to empty string or handle as per your authentication mechanism
                #     administrator=False,
                #     is_admin=False
                # )
                new_user = User()
                new_user.email = email
                new_user.auth_id = auth_id or ''
                new_user.first_name = first_name or ''
                new_user.last_name = last_name or ''
                new_user.institution = institution
                new_user.confirmed_at = datetime.now()
                new_user.is_active = True
                new_user.password = ''  # Set to empty string or handle as per your authentication mechanism
                new_user.administrator = False
                new_user.is_admin = False
                db.session.add(new_user)
                db.session.commit()
                return new_user

        except Exception as e:
            logger.exception(f"Error in UserService: Failed to add or update user: {e}")
            db.session.rollback()  # Rollback session on error
            raise


class UserProjectService:
    @staticmethod
    def add_or_update_user_project(user_id: int, project_id: int, is_owner: bool = False, can_write: bool = False):
        """
        Adds a new user-project relationship or updates an existing one.
        Permission logic:
        - If is_owner is True, can_read and can_write are set to True.
        - If can_write is True (and is_owner is False), can_read is set to True.
        - If neither is_owner nor can_write is True, can_read is set to True.
        """
        try:
            user_project = UserProject.query.filter_by(user_id=user_id, project_id=project_id).first()

            if user_project:
                user_project.is_owner = is_owner

                if is_owner:
                    user_project.can_read = True
                    user_project.can_write = True
                else:
                    user_project.can_write = can_write
                    user_project.can_read = True  # Default to True

                db.session.commit()
                return user_project
            else:
                # Apply permission logic before creating
                if is_owner:
                    can_read = True
                    can_write = True
                elif can_write:
                    can_read = True
                else:
                    can_read = True  # Default to True

                new_user_project = UserProject(
                    user_id=user_id,
                    project_id=project_id,
                    is_owner=is_owner,
                    can_read=can_read,
                    can_write=can_write
                )
                db.session.add(new_user_project)
                db.session.commit()
                return new_user_project

        except Exception as e:
            logger.exception(f"Error in UserProjectService: Failed to add/update user-project entry: {e}")
            db.session.rollback()  # Rollback session on error
            raise

    @staticmethod
    def get_user_project_permissions(user_id: int, project_id: int) -> dict:
        """
        Returns the permission info (can_read, can_write, is_owner) for the given user and project.
        """
        try:
            user_project = UserProject.query.filter_by(user_id=user_id, project_id=project_id).first()
            if not user_project:
                return {"can_read": False, "can_write": False, "is_owner": False}

            return {
                "can_read": user_project.can_read,
                "can_write": user_project.can_write,
                "is_owner": user_project.is_owner
            }

        except Exception as e:
            logger.exception(f"Error in UserProjectService: Failed to get permissions: {e}")
            raise


    @staticmethod
    def remove_user_from_project(user_id: int, project_id: int):
        """Remove a user from a project."""
        try:
            # Fetch the UserProject record
            user_project = UserProject.query.filter_by(user_id=user_id, project_id=project_id).first()
            if not user_project:
                return None  # User is not part of the project
            
            # Delete the record from the database
            db.session.delete(user_project)
            db.session.commit()
            
            logger.info(f"User {user_id} removed from project {project_id}")

        except Exception as e:
            logger.exception(f"Error in remove_user_from_project: {e}")
            db.session.rollback()  # Rollback in case of error
            raise