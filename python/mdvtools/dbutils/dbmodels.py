from flask_sqlalchemy import SQLAlchemy
from datetime import datetime
from mdvtools.dbutils.app import app
from sqlalchemy.orm import relationship

db = SQLAlchemy(app)

class Project(db.Model):
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    name = db.Column(db.String(255), nullable=False, unique=False, default='unnamed_project')
    path = db.Column(db.String(1024), nullable=False, unique=True)
    created_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now)
    is_deleted = db.Column(db.Boolean, nullable=False, default=False)
    deleted_timestamp = db.Column(db.DateTime, nullable=True, default=None)
    update_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    accessed_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now)
    
    def soft_delete(self):
        self.deleted = True
        self.deleted_timestamp = datetime.now()
        db.session.commit()
    
    @classmethod
    def create_project(cls, path):
        new_project = cls(
            path=path,
            created_timestamp=datetime.now(),
            update_timestamp=datetime.now(),
            accessed_timestamp=datetime.now()
        )
        db.session.add(new_project)
        db.session.commit()
        return new_project

class File(db.Model):
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    name = db.Column(db.String(255), nullable=False)
    file_path = db.Column(db.String(255), nullable=True)
    upload_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now)
    update_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    project_id = db.Column(db.Integer, db.ForeignKey('project.id'), nullable=False)

class User(db.Model):
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    username = db.Column(db.String(80), unique=True, nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=False)
    projects = relationship('UserProject', backref='user', lazy=True)

class UserProject(db.Model):
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)
    project_id = db.Column(db.Integer, db.ForeignKey('project.id'), nullable=False)
    can_read = db.Column(db.Boolean, nullable=False, default=False)
    can_write = db.Column(db.Boolean, nullable=False, default=False)

# Function to create default entries
# def create_default_projects():
#     default_projects = ['pbmc3k', 'pbmc3k_project2']
#     with app.app_context():
#         # Create tables if they don't exist
#         db.create_all()

#         # Iterate through default projects
#         for project_name in default_projects:
#             # Check if the project already exists
#             existing_project = Project.query.filter_by(name=project_name).first()

#             # If the project doesn't exist, add it
#             if not existing_project:
#                 project = Project(name=project_name)  # type: ignore
#                 db.session.add(project)

#         # Commit the changes
#         db.session.commit()

# # Call the function to create default entries when the application starts
# create_default_projects()
