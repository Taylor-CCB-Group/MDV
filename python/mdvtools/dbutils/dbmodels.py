from flask_sqlalchemy import SQLAlchemy
from datetime import datetime
from sqlalchemy import String, Integer, DateTime, Boolean, Text, JSON
from sqlalchemy.orm import Mapped, mapped_column, DeclarativeBase, MappedAsDataclass

# https://flask-sqlalchemy.readthedocs.io/en/stable/models/
class Base(DeclarativeBase):
    pass

db = SQLAlchemy(model_class=Base)

class User(db.Model):
    __tablename__ = 'users'
    id = mapped_column(Integer, primary_key=True, autoincrement=True)
    email = mapped_column(String(255), unique=True, nullable=False, default='')
    confirmed_at = mapped_column(DateTime, nullable=True)
    password = mapped_column(String(255), nullable=False, default='')
    is_active = mapped_column(Boolean, nullable=False, default=False)
    first_name = mapped_column(String(50), nullable=False, default='')
    last_name = mapped_column(String(50), nullable=False, default='')
    administrator = mapped_column(Boolean, nullable=False, default=False)
    institution = mapped_column(Text, nullable=True)
    projects = db.relationship('UserProject', backref='user', lazy=True)
    jobs = db.relationship('Job', backref='user', lazy=True)
    permissions = db.relationship('Permission', backref='user', lazy=True)
    preferences = db.relationship('UserPreference', backref='user', lazy=True)
    #shared_objects = db.relationship('SharedObject', foreign_keys='SharedObject.shared_with', backref='shared_with_user', lazy=True)

class Project(db.Model):
    __tablename__ = 'projects'
    id = mapped_column(Integer, primary_key=True, autoincrement=True)
    name = mapped_column(String(255), nullable=False, default='unnamed_project')
    path = mapped_column(String(1024), nullable=False, unique=True)
    created_timestamp = mapped_column(DateTime, nullable=False, default=datetime.now)
    is_deleted = mapped_column(Boolean, nullable=False, default=False)
    deleted_timestamp = mapped_column(DateTime, nullable=True)
    update_timestamp = mapped_column(DateTime, nullable=False, default=datetime.now)
    accessed_timestamp = mapped_column(DateTime, nullable=False, default=datetime.now)
    owner = mapped_column(Integer)
    type = mapped_column(Text)
    data = mapped_column(JSON)
    is_public = mapped_column(Boolean, nullable=False, default=False)
    date_made_public = mapped_column(DateTime)
    status = mapped_column(Text)
    genome = mapped_column(String, db.ForeignKey('genomes.name'))
    parent = mapped_column(Integer)
    description = mapped_column(Text)
    access_level = mapped_column(String(50), nullable=False, default='editable')  # Default access level
    users = db.relationship('UserProject', backref='project', lazy=True)
    files = db.relationship('File', backref='project', lazy=True)

    __table_args__ = (
        db.Index('idx_projects_genome', 'genome'),
        db.Index('idx_projects_owner', 'owner'),
    )


class File(db.Model):
    __tablename__ = 'files'
    id = mapped_column(Integer, primary_key=True, autoincrement=True)
    name = mapped_column(String(255), nullable=False)
    file_path = mapped_column(String(255), nullable=False, unique=True)
    upload_timestamp = mapped_column(DateTime, nullable=False, default=datetime.now)
    update_timestamp = mapped_column(DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    project_id = mapped_column(Integer, db.ForeignKey('projects.id'), nullable=False)
    #project = db.relationship('Project', backref=db.backref('files', lazy=True))
    
class UserProject(db.Model):
    __tablename__ = 'user_projects'
    id = mapped_column(Integer, primary_key=True, autoincrement=True)
    user_id = mapped_column(Integer, db.ForeignKey('users.id'), nullable=False)
    project_id = mapped_column(Integer, db.ForeignKey('projects.id'), nullable=False)
    can_read = mapped_column(Boolean, nullable=False, default=False)
    can_write = mapped_column(Boolean, nullable=False, default=False)

class Genome(db.Model):
    __tablename__ = 'genomes'
    id = mapped_column(Integer, primary_key=True, autoincrement=True)
    name = mapped_column(String(50), nullable=False, unique=True)
    label = mapped_column(Text)
    data = mapped_column(JSON)
    database = mapped_column(Text)
    date_added = mapped_column(DateTime, nullable=False, default=datetime.now)
    connections = mapped_column(Integer)
    icon = mapped_column(Text)
    is_public = mapped_column(Boolean, nullable=False, default=True)
    chrom_sizes = mapped_column(JSON)
    small_icon = mapped_column(Text)
    projects = db.relationship('Project', backref='genomes', lazy=True)

class Job(db.Model):
    __tablename__ = 'jobs'
    id = mapped_column(Integer, primary_key=True, autoincrement=True)
    inputs = mapped_column(JSON)
    user_id = mapped_column(Integer, db.ForeignKey('users.id'), nullable=True)
    outputs = mapped_column(JSON)
    sent_on = mapped_column(DateTime, nullable=False, default=datetime.now)
    status = mapped_column(String(200))
    class_name = mapped_column(String(200))
    genome = mapped_column(String(100))
    finished_on = mapped_column(DateTime)
    is_deleted = mapped_column(Boolean, nullable=False, default=False)
    type = mapped_column(String(200))

class Permission(db.Model):
    __tablename__ = 'permissions'
    id = mapped_column(Integer, primary_key=True, autoincrement=True)
    user_id = mapped_column(Integer, db.ForeignKey('users.id'), nullable=True)
    permission = mapped_column(String(200), nullable=False)
    value = mapped_column(String(200), nullable=False)

class SharedObject(db.Model):
    __tablename__ = 'shared_objects'
    id = mapped_column(Integer, primary_key=True, autoincrement=True)
    owner = mapped_column(Integer)
    shared_with = mapped_column(Integer)
    object_id = mapped_column(Integer)
    date_shared = mapped_column(DateTime, nullable=False, default=datetime.now)
    level = mapped_column(Text, nullable=False, default='view')

class UserPreference(db.Model):
    __tablename__ = 'user_preferences'
    id = mapped_column(Integer, primary_key=True, autoincrement=True)
    preference = mapped_column(Text, nullable=False)
    data = mapped_column(JSON)
    user_id = mapped_column(Integer, db.ForeignKey('users.id'), nullable=True)

class ViewSet(db.Model):
    __tablename__ = 'view_sets'
    id = mapped_column(Integer, primary_key=True, autoincrement=True)
    table_name = mapped_column(String(100))
    name = mapped_column(String(200))
    description = mapped_column(Text)
    date_added = mapped_column(DateTime, nullable=False, default=datetime.now)
    date_modified = mapped_column(DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    owner = mapped_column(Integer, default=0)
    is_public = mapped_column(Boolean, default=False)
    fields = mapped_column(JSON)
    data = mapped_column(JSON)
    date_made_public = mapped_column(DateTime)
    status = mapped_column(Text)
    is_deleted = mapped_column(Boolean, default=False)

class GeneSet(db.Model):
    __tablename__ = 'gene_sets'
    id = mapped_column(Integer, primary_key=True, autoincrement=True)
    name = mapped_column(String(100), nullable=False)
    table_name = mapped_column(String(150))
    data = mapped_column(JSON)
    date_added = mapped_column(DateTime, nullable=False, default=datetime.now)
    date_modified = mapped_column(DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    is_deleted = mapped_column(Boolean, default=False)
    description = mapped_column(Text)

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
