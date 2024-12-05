from flask_sqlalchemy import SQLAlchemy
from datetime import datetime

db = SQLAlchemy()

class User(db.Model):
    __tablename__ = 'users'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    email = db.Column(db.String(255), unique=True, nullable=False, default='')
    confirmed_at = db.Column(db.DateTime, nullable=True)
    password = db.Column(db.String(255), nullable=False, default='')
    is_active = db.Column(db.Boolean, nullable=False, default=False)
    first_name = db.Column(db.String(50), nullable=False, default='')
    last_name = db.Column(db.String(50), nullable=False, default='')
    administrator = db.Column(db.Boolean, nullable=False, default=False)
    institution = db.Column(db.Text, nullable=True)
    projects = db.relationship('UserProject', backref='user', lazy=True)
    jobs = db.relationship('Job', backref='user', lazy=True)
    permissions = db.relationship('Permission', backref='user', lazy=True)
    preferences = db.relationship('UserPreference', backref='user', lazy=True)
    #shared_objects = db.relationship('SharedObject', foreign_keys='SharedObject.shared_with', backref='shared_with_user', lazy=True)

class Project(db.Model):
    __tablename__ = 'projects'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    name = db.Column(db.String(255), nullable=False, default='unnamed_project')
    path = db.Column(db.String(1024), nullable=False, unique=True)
    created_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now)
    is_deleted = db.Column(db.Boolean, nullable=False, default=False)
    deleted_timestamp = db.Column(db.DateTime, nullable=True)
    update_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now)
    accessed_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now)
    owner = db.Column(db.Integer)
    type = db.Column(db.Text)
    data = db.Column(db.JSON)
    is_public = db.Column(db.Boolean, nullable=False, default=False)
    date_made_public = db.Column(db.DateTime)
    status = db.Column(db.Text)
    genome = db.Column(db.String, db.ForeignKey('genomes.name'))
    parent = db.Column(db.Integer)
    description = db.Column(db.Text)
    access_level = db.Column(db.String(50), nullable=False, default='editable')  # Default access level
    users = db.relationship('UserProject', backref='project', lazy=True)
    files = db.relationship('File', backref='project', lazy=True)

    __table_args__ = (
        db.Index('idx_projects_genome', 'genome'),
        db.Index('idx_projects_owner', 'owner'),
    )


class File(db.Model):
    __tablename__ = 'files'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    name = db.Column(db.String(255), nullable=False)
    file_path = db.Column(db.String(255), nullable=False, unique=True)
    upload_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now)
    update_timestamp = db.Column(db.DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    project_id = db.Column(db.Integer, db.ForeignKey('projects.id'), nullable=False)
    #project = db.relationship('Project', backref=db.backref('files', lazy=True))
    
class UserProject(db.Model):
    __tablename__ = 'user_projects'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'), nullable=False)
    project_id = db.Column(db.Integer, db.ForeignKey('projects.id'), nullable=False)
    can_read = db.Column(db.Boolean, nullable=False, default=False)
    can_write = db.Column(db.Boolean, nullable=False, default=False)

class Genome(db.Model):
    __tablename__ = 'genomes'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    name = db.Column(db.String(50), nullable=False, unique=True)
    label = db.Column(db.Text)
    data = db.Column(db.JSON)
    database = db.Column(db.Text)
    date_added = db.Column(db.DateTime, nullable=False, default=datetime.now)
    connections = db.Column(db.Integer)
    icon = db.Column(db.Text)
    is_public = db.Column(db.Boolean, nullable=False, default=True)
    chrom_sizes = db.Column(db.JSON)
    small_icon = db.Column(db.Text)
    projects = db.relationship('Project', backref='genomes', lazy=True)

class Job(db.Model):
    __tablename__ = 'jobs'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    inputs = db.Column(db.JSON)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'), nullable=True)
    outputs = db.Column(db.JSON)
    sent_on = db.Column(db.DateTime, nullable=False, default=datetime.now)
    status = db.Column(db.String(200))
    class_name = db.Column(db.String(200))
    genome = db.Column(db.String(100))
    finished_on = db.Column(db.DateTime)
    is_deleted = db.Column(db.Boolean, nullable=False, default=False)
    type = db.Column(db.String(200))

class Permission(db.Model):
    __tablename__ = 'permissions'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'), nullable=True)
    permission = db.Column(db.String(200), nullable=False)
    value = db.Column(db.String(200), nullable=False)

class SharedObject(db.Model):
    __tablename__ = 'shared_objects'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    owner = db.Column(db.Integer)
    shared_with = db.Column(db.Integer)
    object_id = db.Column(db.Integer)
    date_shared = db.Column(db.DateTime, nullable=False, default=datetime.now)
    level = db.Column(db.Text, nullable=False, default='view')

class UserPreference(db.Model):
    __tablename__ = 'user_preferences'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    preference = db.Column(db.Text, nullable=False)
    data = db.Column(db.JSON)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'), nullable=True)

class ViewSet(db.Model):
    __tablename__ = 'view_sets'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    table_name = db.Column(db.String(100))
    name = db.Column(db.String(200))
    description = db.Column(db.Text)
    date_added = db.Column(db.DateTime, nullable=False, default=datetime.now)
    date_modified = db.Column(db.DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    owner = db.Column(db.Integer, default=0)
    is_public = db.Column(db.Boolean, default=False)
    fields = db.Column(db.JSON)
    data = db.Column(db.JSON)
    date_made_public = db.Column(db.DateTime)
    status = db.Column(db.Text)
    is_deleted = db.Column(db.Boolean, default=False)

class GeneSet(db.Model):
    __tablename__ = 'gene_sets'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True)
    name = db.Column(db.String(100), nullable=False)
    table_name = db.Column(db.String(150))
    data = db.Column(db.JSON)
    date_added = db.Column(db.DateTime, nullable=False, default=datetime.now)
    date_modified = db.Column(db.DateTime, nullable=False, default=datetime.now, onupdate=datetime.now)
    is_deleted = db.Column(db.Boolean, default=False)
    description = db.Column(db.Text)

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
