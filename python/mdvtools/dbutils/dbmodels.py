from flask_sqlalchemy import SQLAlchemy
from datetime import datetime
from sqlalchemy import String, Integer, DateTime, Boolean, Text, JSON, ForeignKey
from sqlalchemy.orm import Mapped, mapped_column, DeclarativeBase, relationship
from typing import Any
# https://flask-sqlalchemy.readthedocs.io/en/stable/models/
class Base(DeclarativeBase):
    pass

db = SQLAlchemy(model_class=Base)

class User(db.Model):
    __tablename__ = 'users'
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    email: Mapped[str] = mapped_column(String(255), unique=True, nullable=False, default='')
    confirmed_at: Mapped[datetime] = mapped_column(DateTime, nullable=True)
    password: Mapped[str] = mapped_column(String(255), nullable=False, default='')
    is_active: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    first_name: Mapped[str] = mapped_column(String(50), nullable=False, default='')
    last_name: Mapped[str] = mapped_column(String(50), nullable=False, default='')
    administrator: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    institution: Mapped[str] = mapped_column(Text, nullable=True)
    projects = relationship('UserProject', backref='user', lazy=True)
    jobs = relationship('Job', backref='user', lazy=True)
    permissions = relationship('Permission', backref='user', lazy=True)
    preferences = relationship('UserPreference', backref='user', lazy=True)
    #shared_objects = relationship('SharedObject', foreign_keys='SharedObject.shared_with', backref='shared_with_user', lazy=True)

class Project(db.Model):
    __tablename__ = 'projects'
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False, default='unnamed_project')
    path: Mapped[str] = mapped_column(String(1024), nullable=False, unique=True)
    created_timestamp: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now)
    is_deleted: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    deleted_timestamp: Mapped[datetime] = mapped_column(DateTime, nullable=True)
    update_timestamp: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now)
    accessed_timestamp: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now)
    owner: Mapped[int] = mapped_column(Integer)
    type: Mapped[str] = mapped_column(Text)
    data: Mapped[Any] = mapped_column(JSON)
    is_public: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    date_made_public: Mapped[datetime] = mapped_column(DateTime)
    status: Mapped[str] = mapped_column(Text)
    genome: Mapped[str] = mapped_column(String, ForeignKey('genomes.name'))
    parent: Mapped[str] = mapped_column(Integer)
    description: Mapped[str] = mapped_column(Text)
    access_level: Mapped[str] = mapped_column(String(50), nullable=False, default='editable')  # Default access level
    users = relationship('UserProject', backref='project', lazy=True)
    files = relationship('File', backref='project', lazy=True)

    __table_args__ = (
        db.Index('idx_projects_genome', 'genome'),
        db.Index('idx_projects_owner', 'owner'),
    )


class File(db.Model):
    __tablename__ = 'files'
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    file_path: Mapped[str] = mapped_column(String(255), nullable=False, unique=True)
    upload_timestamp: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now)
    update_timestamp: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    project_id: Mapped[int] = mapped_column(Integer, ForeignKey('projects.id'), nullable=False)
    #project = relationship('Project', backref=db.backref('files', lazy=True))
    
class UserProject(db.Model):
    __tablename__ = 'user_projects'
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    user_id: Mapped[int] = mapped_column(Integer, ForeignKey('users.id'), nullable=False)
    project_id: Mapped[int] = mapped_column(Integer, ForeignKey('projects.id'), nullable=False)
    can_read: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    can_write: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)

class Genome(db.Model):
    __tablename__ = 'genomes'
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(50), nullable=False, unique=True)
    label: Mapped[str] = mapped_column(Text)
    data: Mapped[Any] = mapped_column(JSON)
    database: Mapped[str] = mapped_column(Text)
    date_added: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now)
    connections: Mapped[int] = mapped_column(Integer)
    icon: Mapped[str] = mapped_column(Text)
    is_public: Mapped[bool] = mapped_column(Boolean, nullable=False, default=True)
    chrom_sizes: Mapped[Any] = mapped_column(JSON)
    small_icon: Mapped[str] = mapped_column(Text)
    projects = relationship('Project', backref='genomes', lazy=True)

class Job(db.Model):
    __tablename__ = 'jobs'
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    inputs: Mapped[Any] = mapped_column(JSON)
    user_id: Mapped[int] = mapped_column(Integer, ForeignKey('users.id'), nullable=True)
    outputs: Mapped[Any] = mapped_column(JSON)
    sent_on: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now)
    status: Mapped[str] = mapped_column(String(200))
    class_name: Mapped[str] = mapped_column(String(200))
    genome: Mapped[str] = mapped_column(String(100))
    finished_on: Mapped[datetime] = mapped_column(DateTime)
    is_deleted: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    type: Mapped[str] = mapped_column(String(200))

class Permission(db.Model):
    __tablename__ = 'permissions'
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    user_id: Mapped[int] = mapped_column(Integer, ForeignKey('users.id'), nullable=True)
    permission: Mapped[str] = mapped_column(String(200), nullable=False)
    value: Mapped[str] = mapped_column(String(200), nullable=False)

class SharedObject(db.Model):
    __tablename__ = 'shared_objects'
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    owner: Mapped[int] = mapped_column(Integer)
    shared_with: Mapped[int] = mapped_column(Integer)
    object_id: Mapped[int] = mapped_column(Integer)
    date_shared: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now)
    level: Mapped[str] = mapped_column(Text, nullable=False, default='view')

class UserPreference(db.Model):
    __tablename__ = 'user_preferences'
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    preference: Mapped[str] = mapped_column(Text, nullable=False)
    data: Mapped[Any] = mapped_column(JSON)
    user_id: Mapped[int] = mapped_column(Integer, ForeignKey('users.id'), nullable=True)

class ViewSet(db.Model):
    __tablename__ = 'view_sets'
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    table_name: Mapped[str] = mapped_column(String(100))
    name: Mapped[str] = mapped_column(String(200))
    description: Mapped[str] = mapped_column(Text)
    date_added: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now)
    date_modified: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    owner: Mapped[int] = mapped_column(Integer, default=0)
    is_public: Mapped[bool] = mapped_column(Boolean, default=False)
    fields: Mapped[Any] = mapped_column(JSON)
    data: Mapped[Any] = mapped_column(JSON)
    date_made_public: Mapped[datetime] = mapped_column(DateTime)
    status: Mapped[str] = mapped_column(Text)
    is_deleted: Mapped[bool] = mapped_column(Boolean, default=False)

class GeneSet(db.Model):
    __tablename__ = 'gene_sets'
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(100), nullable=False)
    table_name: Mapped[str] = mapped_column(String(150))
    data: Mapped[Any] = mapped_column(JSON)
    date_added: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now)
    date_modified: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    is_deleted: Mapped[bool] = mapped_column(Boolean, default=False)
    description: Mapped[str] = mapped_column(Text)

# Indexes
db.Index('idx_genes_name', GeneSet.name)
db.Index('idx_views_table_name', ViewSet.table_name)
db.Index('idx_views_name', ViewSet.name)

"""
@db.event.listens_for(File, 'after_insert')
@db.event.listens_for(File, 'after_update')
@db.event.listens_for(File, 'after_delete')
def update_project_timestamp(mapper, connection, target):
    try:
        # Ensure target.project is the correct way to access the associated Project
        print("--------Event Listener on File")
        project = target.project
        print(project)
        if project:
            # Update the project's timestamp
            project.update_timestamp = datetime.now()
            db.session.commit()
            print(f"Updated project timestamp for project ID {project.id}")
        else:
            print(f"No associated project found for file ID {target.id}")
    except Exception as e:
        print(f"Error updating project timestamp: {e}")
        db.session.rollback()  # Rollback in case of an error

"""
