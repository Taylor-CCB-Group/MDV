"""
Database Service Layer for MDV

This module provides business logic services for database operations related to
projects, files, users, and user-project relationships. All services use static
methods and handle transaction management with proper error handling and rollback.

Classes:
    ProjectService: Manages project-related database operations
    FileService: Manages file-related database operations
    UserService: Manages user-related database operations
    UserProjectService: Manages user-project relationships and permissions

All services follow a consistent pattern:
    - Static methods for stateless operations
    - Automatic transaction management (commit on success, rollback on error)
    - Comprehensive error logging
    - Proper exception propagation

Example:
    >>> from mdvtools.dbutils.dbservice import ProjectService
    >>> projects = ProjectService.get_active_projects()
    >>> new_project = ProjectService.add_new_project('/path/to/project', 'My Project')
"""

from mdvtools.dbutils.dbmodels import db, Project, File, User, UserProject
from datetime import datetime
from mdvtools.mdvproject import MDVProject
from typing import Optional
from mdvtools.logging_config import get_logger

logger = get_logger(__name__)


class ProjectService:
    """
    Service class for managing project-related database operations.
    
    This service provides methods for creating, retrieving, updating, and
    deleting projects in the database. It handles project metadata, timestamps,
    and access levels.
    
    Class Attributes:
        failed_projects (list[tuple[int, str | Exception]]): List of tuples
            containing failed project IDs and associated error messages/exceptions.
            Projects in this list are excluded from get_active_projects() results.
    
    All methods are static and handle their own transaction management.
    """
    
    # list of tuples containing failed project IDs and associated error messages/exceptions
    failed_projects: list[tuple[int, str | Exception]] = []
    
    @staticmethod
    def get_active_projects():
        """
        Retrieve all active (non-deleted) projects with their metadata.
        
        Returns a list of dictionaries containing project information including
        ID, name, last modified timestamp, thumbnail image, and README content.
        Projects that are in the failed_projects list are excluded from results.
        
        Returns:
            list[dict]: List of project dictionaries, each containing:
                - id (int): Project ID
                - name (str): Project name
                - lastModified (str): Formatted update timestamp (YYYY-MM-DD HH:MM:SS)
                - thumbnail (str|None): First available viewImage from project views
                - readme (str|None): Project README content if available
        
        Raises:
            Exception: If database query fails or project loading encounters errors.
            The exception is logged and re-raised.
        
        Note:
            This method is used by the /projects API route.
            Thumbnail extraction may fail silently and return None if the project
            cannot be loaded or has no viewImage.
        """
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
    def get_next_project_id():
        """
        Calculate the next available project ID.
        
        Queries the database for the maximum existing project ID and returns
        the next sequential ID. If no projects exist, returns 1.
        
        Returns:
            int: The next available project ID (minimum 1)
        
        Raises:
            Exception: If database query fails. The exception is logged and re-raised.
        
        Note:
            This method is used by the /create_project API route.
        """
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
    def add_new_project(path, name='unnamed_project'):
        """
        Create a new project record in the database.
        
        Creates a new Project instance with the specified path and name,
        adds it to the database session, and commits the transaction.
        
        Args:
            path (str): Filesystem path to the project directory. Must be unique.
            name (str, optional): Project name. Defaults to 'unnamed_project'.
        
        Returns:
            Project: The newly created Project model instance.
        
        Raises:
            Exception: If database operation fails (e.g., duplicate path).
            The exception is logged, transaction is rolled back, and exception
            is re-raised.
        
        Note:
            This method is used by the /create_project API route.
            The project's created_timestamp, update_timestamp, and
            accessed_timestamp are automatically set to the current datetime.
        """
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
    def get_project_by_id(id):
        """
        Retrieve a project by its ID.
        
        Args:
            id (int): The project ID to look up.
        
        Returns:
            Project|None: The Project model instance if found, None otherwise.
        
        Raises:
            Exception: If database query fails. The exception is logged and re-raised.
        
        Note:
            This method is used by the /delete_project API route and other
            internal methods that need to fetch a project by ID.
        """
        try:
            return Project.query.get(id)
        except Exception as e:
            logger.exception(f"Error in dbservice: Error querying project by id-No project found with id: {e}")
            raise

    @staticmethod
    def soft_delete_project(id):
        """
        Perform a soft delete on a project.
        
        Sets the project's is_deleted flag to True and records the deletion
        timestamp. The project record remains in the database but is excluded
        from normal queries.
        
        Args:
            id (int): The project ID to soft delete.
        
        Returns:
            bool: True if the project was found and deleted, False if the
                project does not exist.
        
        Raises:
            Exception: If database operation fails. The exception is logged,
            transaction is rolled back, and exception is re-raised.
        
        Note:
            This method is used by the /delete_project API route.
            Soft deletion preserves data integrity and allows for potential
            recovery or audit purposes.
        """
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
    def update_project_name(project_id, new_name):
        """
        Update a project's name.
        
        Changes the project name and updates both the update_timestamp and
        accessed_timestamp to the current datetime.
        
        Args:
            project_id (int): The project ID to update.
            new_name (str): The new name for the project.
        
        Returns:
            bool: True if the project was found and updated, False if the
                project does not exist.
        
        Raises:
            Exception: If database operation fails. The exception is logged,
            transaction is rolled back, and exception is re-raised.
        
        Note:
            This method is used by the /rename_project API route.
        """
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
    def change_project_access(project_id, new_access_level):
        """
        Change the access level of a project.
        
        Updates the project's access_level field and updates both the
        update_timestamp and accessed_timestamp to the current datetime.
        
        Args:
            project_id (int): The project ID to update.
            new_access_level (str): The new access level (e.g., 'editable',
                'readonly', 'private').
        
        Returns:
            tuple: A tuple containing:
                - str|None: The new access level if successful, None if project not found
                - str: Status message ('Success' or error message)
                - int: HTTP status code (200 for success, 404 if project not found)
        
        Raises:
            Exception: If database operation fails. The exception is logged,
            transaction is rolled back, and exception is re-raised.
        
        Note:
            This method is used by the /access API route.
        """
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
    def set_project_update_timestamp(project_id: str):
        """
        Set the project's update timestamp to the current datetime.
        
        Updates the update_timestamp field for a project, typically called
        when project data or metadata is modified.
        
        Args:
            project_id (str|int): The project ID to update. Can be string or int.
        
        Returns:
            None: Returns silently if project is not found (logged as info).
        
        Raises:
            Exception: If database operation fails. The exception is logged,
            transaction is rolled back, and exception is re-raised.
        
        Note:
            This method is typically called internally when project data changes.
            It does not update the accessed_timestamp.
        """
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
    def set_project_accessed_timestamp(project_id: str):
        """
        Set the project's accessed timestamp to the current datetime.
        
        Updates the accessed_timestamp field for a project, typically called
        when a project is viewed or accessed by a user.
        
        Args:
            project_id (str|int): The project ID to update. Can be string or int.
        
        Returns:
            None: Returns silently if project is not found (logged as info).
        
        Raises:
            Exception: If database operation fails. The exception is logged,
            transaction is rolled back, and exception is re-raised.
        
        Note:
            This method is typically called when a project is accessed for
            tracking purposes. It does not update the update_timestamp.
        """
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
    """
    Service class for managing file-related database operations.
    
    This service provides methods for adding, updating, retrieving, and
    deleting file records associated with projects. Files are tracked with
    their paths, names, and timestamps.
    
    All methods are static and handle their own transaction management.
    """
    
    @staticmethod
    def add_or_update_file_in_project(file_name, file_path, project_id):
        """
        Add a new file or update an existing file in the database.
        
        If a file with the same file_path and project_id already exists,
        updates its name if different. Otherwise, creates a new File record.
        
        Args:
            file_name (str): The name of the file.
            file_path (str): The filesystem path to the file. Must be unique
                per project.
            project_id (int): The ID of the project this file belongs to.
        
        Returns:
            None: This method commits the transaction but does not return a value.
        
        Raises:
            Exception: If database operation fails. The exception is logged,
            transaction is rolled back, and exception is re-raised.
        
        Note:
            New files have both upload_timestamp and update_timestamp set to
            the current datetime. Updated files have only update_timestamp updated.
        """
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
        """
        Fetch a file by its path and project ID.
        
        Args:
            file_path (str): The filesystem path to the file.
            project_id (int): The ID of the project this file belongs to.
        
        Returns:
            File|None: The File model instance if found, None if not found
                or if an error occurs.
        
        Note:
            Errors are logged but do not raise exceptions. Returns None on error.
        """
        try:
            return File.query.filter_by(file_path=file_path, project_id=project_id).first()
        except Exception as e:
            logger.exception(f"Error retrieving file by path '{file_path}' and project ID {project_id}: {e}")
            return None
    
    @staticmethod
    def file_exists_in_project(file_path, project_id):
        """
        Check if a file exists in the database for a given project.
        
        Utility function to determine if a file record exists without
        retrieving the full File object.
        
        Args:
            file_path (str): The filesystem path to the file.
            project_id (int): The ID of the project this file belongs to.
        
        Returns:
            bool: True if the file exists, False if it doesn't exist or if
                an error occurs.
        
        Note:
            Errors are logged but do not raise exceptions. Returns False on error.
        """
        try:
            return File.query.filter_by(
                file_path=file_path, project_id=project_id
            ).first() is not None
        except Exception as e:
            logger.exception(f"Error checking file existence: {e}")
            return False

    @staticmethod
    def get_files_by_project(project_id):
        """
        Retrieve all files associated with a project.
        
        Args:
            project_id (int): The ID of the project to get files for.
        
        Returns:
            list[File]: List of File model instances for the project.
                Returns empty list if no files found or if an error occurs.
        
        Note:
            Errors are logged but do not raise exceptions. Returns empty list on error.
        """
        try:
            return File.query.filter_by(project_id=project_id).all()
        except Exception as e:
            logger.exception(f"Error querying files for project ID {project_id}: {e}")
            return []

    @staticmethod
    def delete_files_by_project(project_id):
        """
        Delete all files associated with a project.
        
        Permanently removes all File records for the specified project from
        the database. This is typically used when deleting a project.
        
        Args:
            project_id (int): The ID of the project whose files should be deleted.
        
        Returns:
            bool: True if the operation completed successfully, False if an
                error occurred.
        
        Note:
            Errors are logged and transaction is rolled back. Returns False on error.
            This is a permanent deletion, not a soft delete.
        """
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
        """
        Update the update_timestamp for a file.
        
        Sets the file's update_timestamp to the current datetime. Typically
        called when file metadata or content is modified.
        
        Args:
            file_id (int): The ID of the file to update.
        
        Returns:
            bool: True if the file was found and updated, False if the file
                does not exist or if an error occurred.
        
        Note:
            Errors are logged and transaction is rolled back. Returns False on error.
        """
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
    """
    Service class for managing user-related database operations.
    
    This service provides methods for creating and updating user records.
    Users are identified by their email address, which must be unique.
    
    All methods are static and handle their own transaction management.
    """
    
    @staticmethod
    def add_or_update_user(email: str, auth_id: Optional[str] = None, first_name: Optional[str] = None, last_name: Optional[str] = None, institution: Optional[str] = None):
        """
        Add a new user or update an existing user based on email address.
        
        If a user with the given email exists, updates the provided fields.
        If no user exists, creates a new user record with the provided information.
        Email is required and must be unique.
        
        Args:
            email (str): User's email address. Required and must be unique.
            auth_id (str, optional): User's authentication provider ID (e.g., Auth0).
            first_name (str, optional): User's first name.
            last_name (str, optional): User's last name.
            institution (str, optional): User's institution or organization.
        
        Returns:
            User: The created or updated User model instance.
        
        Raises:
            ValueError: If email is not provided or is empty.
            Exception: If database operation fails. The exception is logged,
                transaction is rolled back, and exception is re-raised.
        
        Note:
            New users are created with:
            - confirmed_at set to current datetime
            - is_active set to True
            - password set to empty string
            - administrator and is_admin set to False
            - Empty strings for optional fields if not provided
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
    """
    Service class for managing user-project relationships and permissions.
    
    This service provides methods for managing access control between users
    and projects. It handles permission levels including read, write, and
    ownership permissions.
    
    Permission Logic:
        - If is_owner is True: can_read and can_write are automatically True
        - If can_write is True (and is_owner is False): can_read is set to True
        - If neither is_owner nor can_write is True: can_read defaults to True
    
    All methods are static and handle their own transaction management.
    """
    
    @staticmethod
    def add_or_update_user_project(user_id: int, project_id: int, is_owner: bool = False, can_write: bool = False):
        """
        Add a new user-project relationship or update an existing one.
        
        Creates or updates the permission relationship between a user and a project.
        Permission levels are automatically set according to the service's
        permission logic (see class docstring).
        
        Args:
            user_id (int): The ID of the user.
            project_id (int): The ID of the project.
            is_owner (bool, optional): Whether the user is the project owner.
                Defaults to False. If True, automatically grants read and write.
            can_write (bool, optional): Whether the user can write to the project.
                Defaults to False. If True (and is_owner is False), grants read access.
        
        Returns:
            UserProject: The created or updated UserProject model instance.
        
        Raises:
            Exception: If database operation fails. The exception is logged,
                transaction is rolled back, and exception is re-raised.
        
        Note:
            When updating an existing relationship, all permission fields are
            recalculated based on the provided is_owner and can_write values.
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
        Get the permission information for a user-project relationship.
        
        Retrieves the read, write, and ownership permissions for a specific
        user and project combination.
        
        Args:
            user_id (int): The ID of the user.
            project_id (int): The ID of the project.
        
        Returns:
            dict: A dictionary containing permission flags:
                - can_read (bool): Whether the user can read the project
                - can_write (bool): Whether the user can write to the project
                - is_owner (bool): Whether the user owns the project
        
            If no relationship exists, returns all flags as False.
        
        Raises:
            Exception: If database query fails. The exception is logged and re-raised.
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
        """
        Remove a user from a project by deleting the user-project relationship.
        
        Permanently removes the UserProject record, effectively revoking all
        permissions for the user on the specified project.
        
        Args:
            user_id (int): The ID of the user to remove.
            project_id (int): The ID of the project to remove the user from.
        
        Returns:
            None: Returns None if the user was not part of the project (no-op).
                Otherwise commits the deletion.
        
        Raises:
            Exception: If database operation fails. The exception is logged,
                transaction is rolled back, and exception is re-raised.
        
        Note:
            This is a permanent deletion. If the user needs access again,
            add_or_update_user_project() must be called to recreate the relationship.
        """
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