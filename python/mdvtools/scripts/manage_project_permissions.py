#!/usr/bin/env python3
import argparse
import sys
import os
from pathlib import Path
from typing import Dict, Union

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

def batch_assign_from_text(file_path: str) -> bool:
    """
    Batch assign permissions from a plain text file.

    Accepted line formats (comments and empty lines are ignored):
        project_name,email,permission
        project_name email permission
    """
    path = Path(file_path)
    if not path.exists():
        print(f"Error: file not found: {file_path}")
        return False

    with path.open(encoding="utf-8") as f:
        lines = [
            line.strip()
            for line in f
            if line.strip() and not line.strip().startswith("#")
        ]

    if not lines:
        print("No assignments found in file.")
        return False

    allowed_perms = {"view", "edit", "owner"}
    success = True

    for idx, line in enumerate(lines, start=1):
        # Support comma-separated or whitespace-separated
        if "," in line:
            parts = [p.strip() for p in line.split(",")]
        else:
            parts = line.split()

        if len(parts) != 3:
            print(
                f"Line {idx}: invalid format. Expected 'project_name,email,permission' "
                f"or 'project_name email permission', got: {line}"
            )
            success = False
            continue

        project_name, email, permission = parts

        if permission not in allowed_perms:
            print(
                f"Line {idx}: invalid permission '{permission}'. "
                f"Must be one of: view, edit, owner."
            )
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
    
    # Batch assignment command (text only)
    batch_txt_parser = subparsers.add_parser(
        'batch_txt',
        help='Batch assign permissions from a plain text file '
             '(lines: project_name,email,permission or project_name email permission)',
    )
    batch_txt_parser.add_argument('--file', required=True, help='Path to text file with assignments')

    args = parser.parse_args()

    # Create application context
    with app.app_context():
        if args.command == 'assign':
            success = assign_permissions(args.email, args.project, args.permission)
        elif args.command == 'batch_txt':
            success = batch_assign_from_text(args.file)
        else:
            parser.print_help()
            sys.exit(1)
    
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main() 