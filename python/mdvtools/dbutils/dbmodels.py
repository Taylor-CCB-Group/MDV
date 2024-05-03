from mdvtools.dbutils.app import app
from datetime import datetime
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy.orm import relationship
from dataclasses import dataclass

db = SQLAlchemy(app)

@dataclass
class Project(db.Model):
    id: int
    name: str
    #files: 'list[File]'  # Remove the type annotation here

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(255), nullable=False, unique=True)
    files = relationship('File', backref='project', lazy=True)

@dataclass
class File(db.Model):
    id: int
    name: str
    file_path: str
    upload_timestamp: datetime
    update_timestamp: datetime
    project_id: int

    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(255), nullable=False)
    file_path = db.Column(db.String(255), nullable=False)
    upload_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now)
    update_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    project_id = db.Column(db.Integer, db.ForeignKey('project.id'), nullable=False)

@dataclass
class User(db.Model):
    id: int
    username: str
    email: str

    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(80), unique=True, nullable=False)
    email = db.Column(db.String(120), unique=True, nullable=False)
