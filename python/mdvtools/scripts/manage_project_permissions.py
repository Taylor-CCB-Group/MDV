#!/usr/bin/env python3
import argparse
import sys
import os
import json
from typing import List, Dict, Union

# Add the parent directory to sys.path to import mdvtools modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from mdvtools.dbutils.dbmodels import db, User, Project, UserProject
from mdvtools.dbutils.dbservice import UserProjectService
from mdvtools.auth.authutils import cache_user_projects
from mdvtools.dbutils.safe_mdv_app import app  # Import the Flask app

def get_user_by_email(email: str) -> Union[User, None]:
    """Get user by email."""
    return User.query.filter_by(email=email).first()

def get_project_by_name(name: str) -> Union[Project, None]:
    """Get project by name."""
    return Project.query.filter_by(name=name, is_deleted=False).first()

def parse_permission(permission: str) -> Dict[str, bool]:
    """Parse permission string into a dict of boolean flags."""
    return {
        'is_owner': permission == 'owner',
        'can_write': permission in ['edit', 'owner'],
        'can_read': True  # Always true if added to a project
    }

def assign_permissions(user_email: str, project_name: str, permission: str):
    """Assign permissions for a single user and project."""
    user = get_user_by_email(user_email)
    if not user:
        print(f"Error: User with email {user_email} not found")
        return False

    project = get_project_by_name(project_name)
    if not project:
        print(f"Error: Project with name {project_name} not found")
        return False

    perm = parse_permission(permission)
    try:
        # Update database
        UserProjectService.add_or_update_user_project(
            user_id=user.id,
            project_id=project.id,
            is_owner=perm['is_owner'],
            can_write=perm['can_write']
        )

        # Refresh all caches
        cache_user_projects()

        print(f"Successfully assigned {permission} permission to {user_email} for project {project_name}")
        return True
    except Exception as e:
        print(f"Error assigning permissions: {str(e)}")
        return False

def batch_assign_from_file(file_path: str):
    """
    Batch assign permissions from a JSON file.
    Expected format:
    {
        "assignments": [
            {
                "email": "user@example.com",
                "projects": [
                    {"name": "project1", "permission": "view"},
                    {"name": "project2", "permission": "edit"}
                ]
            }
        ]
    }
    """
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading file: {str(e)}")
        return False

    success = True
    for assignment in data.get('assignments', []):
        email = assignment.get('email')
        if not email:
            print("Error: Missing email in assignment")
            success = False
            continue

        for project in assignment.get('projects', []):
            project_name = project.get('name')
            permission = project.get('permission')
            if not project_name or not permission:
                print(f"Error: Missing project name or permission for {email}")
                success = False
                continue

            if not assign_permissions(email, project_name, permission):
                success = False

    return success

def main():
    parser = argparse.ArgumentParser(description='Manage project permissions for users')
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Single assignment command
    assign_parser = subparsers.add_parser('assign', help='Assign permission for a single user and project')
    assign_parser.add_argument('--email', required=True, help='User email')
    assign_parser.add_argument('--project', required=True, help='Project name')
    assign_parser.add_argument('--permission', required=True, choices=['view', 'edit', 'owner'], 
                             help='Permission level')
    
    # Batch assignment command
    batch_parser = subparsers.add_parser('batch', help='Batch assign permissions from a JSON file')
    batch_parser.add_argument('--file', required=True, help='Path to JSON file with assignments')

    args = parser.parse_args()

    # Create application context
    with app.app_context():
        if args.command == 'assign':
            success = assign_permissions(args.email, args.project, args.permission)
        elif args.command == 'batch':
            success = batch_assign_from_file(args.file)
        else:
            parser.print_help()
            sys.exit(1)
    
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main() 