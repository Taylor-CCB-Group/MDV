# project_service.py

from mdvtools.dbutils.dbmodels import db, Project
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
