from dataclasses import dataclass, field
from typing import Optional
from mdvtools.dbutils.app import app
from datetime import datetime
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy.orm import relationship

db = SQLAlchemy(app)

@dataclass
class Project(db.Model):
    id: int = field(init=False)
    name: str
    #files: 'list[File]' = field(default_factory=list)

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    name = db.Column(db.String(255), nullable=False, unique=True)
    files = relationship('File', backref='project', lazy=True)

@dataclass
class File(db.Model):
    id: int = field(init=False)
    name: str
    file_path: Optional[str] = None
    upload_timestamp: datetime = field(default_factory=datetime.now)
    update_timestamp: datetime = field(default_factory=datetime.now, compare=False)
    project_id: int = field(init=False)

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    name = db.Column(db.String(255), nullable=False)
    file_path = db.Column(db.String(255), nullable=True)
    upload_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now)
    update_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    project_id = db.Column(db.Integer, db.ForeignKey('project.id'), nullable=False, autoincrement=True)

@dataclass
class User(db.Model):
    id: int = field(init=False)
    username: str
    email: str

    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    username = db.Column(db.String(80), unique=True, nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=False)


# Function to create default entries
def create_default_projects():
    default_projects = ['pbmc3k', 'pbmc3k_project2']
    with app.app_context():
        # Create tables if they don't exist
        db.create_all()
        
        # Iterate through default projects
        for project_name in default_projects:
            # Check if the project already exists
            existing_project = Project.query.filter_by(name=project_name).first()
            
            # If the project doesn't exist, add it
            if not existing_project:
                project = Project(name=project_name)
                db.session.add(project)
        
        # Commit the changes
        db.session.commit()



# Call the function to create default entries when the application starts
create_default_projects()