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

def sync_users_from_auth0():
    """Sync users from Auth0 to the database if Auth0 is enabled."""
    try:
        from mdvtools.auth.authutils import get_auth_provider
        auth_provider = get_auth_provider()
        if hasattr(auth_provider, 'sync_users_to_db'):
            print("Syncing users from Auth0 to database...")
            auth_provider.sync_users_to_db()
            print("User sync completed.")
            return True
        else:
            print("Auth provider does not support user sync.")
            return False
    except Exception as e:
        print(f"Warning: Could not sync users from Auth0: {e}")
        return False

def assign_permissions(user_email: str, project_name: str, permission: str, refresh_cache: bool = True):
    """
    Assign permissions for a single user and project.
    
    Args:
        user_email: User's email address
        project_name: Project name
        permission: Permission level (view, edit, owner)
        refresh_cache: Whether to refresh cache after assignment (default: True)
    """
    user = get_user_by_email(user_email)
    if not user:
        # Try syncing from Auth0 first
        print(f"User with email {user_email} not found in database. Syncing from Auth0...")
        sync_users_from_auth0()
        if refresh_cache:
            cache_user_projects()
        # Try again after sync
        user = get_user_by_email(user_email)
        if not user:
            print(f"Error: User with email {user_email} not found even after sync")
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

        # Refresh cache only if requested (skip in batch mode for efficiency)
        if refresh_cache:
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

    # Extract all unique emails from the file to check if they exist
    emails_to_check = set()
    allowed_perms = {"view", "edit", "owner"}
    success = True

    # First pass: parse and validate format, collect emails
    parsed_lines = []
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

        emails_to_check.add(email)
        parsed_lines.append((project_name, email, permission))

    # Check if any users are missing and sync from Auth0 if needed
    missing_users = []
    for email in emails_to_check:
        user = get_user_by_email(email)
        if not user:
            missing_users.append(email)

    if missing_users:
        print(f"Found {len(missing_users)} user(s) not in database. Syncing from Auth0...")
        sync_users_from_auth0()
        # Refresh cache after sync
        cache_user_projects()

    # Second pass: assign permissions (without individual cache refreshes for efficiency)
    for project_name, email, permission in parsed_lines:
        if not assign_permissions(email, project_name, permission, refresh_cache=False):
            success = False

    # Final cache refresh to ensure all changes are reflected
    print("Refreshing cache...")
    cache_user_projects()
    print("Cache refreshed. All permission changes are now active.")

    return success

def main():
    # Check if first argument is a file path (before parsing, to avoid subparser conflicts)
    if len(sys.argv) > 1:
        first_arg = sys.argv[1]
        # If it looks like a file path (contains / or .txt or exists as a file), treat it as such
        if '/' in first_arg or first_arg.endswith('.txt') or Path(first_arg).exists():
            # Direct file path mode - process immediately
            with app.app_context():
                success = batch_assign_from_text(first_arg)
            sys.exit(0 if success else 1)
    
    # Otherwise, use normal command parsing for single assignment
    parser = argparse.ArgumentParser(
        description='Manage project permissions for users',
        epilog='Usage: python manage_project_permissions.py <file_path> OR python manage_project_permissions.py assign --email ... --project ... --permission ...'
    )
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', help='Commands', required=True)
    
    # Single assignment command
    assign_parser = subparsers.add_parser('assign', help='Assign permission for a single user and project')
    assign_parser.add_argument('--email', required=True, help='User email')
    assign_parser.add_argument('--project', required=True, help='Project name')
    assign_parser.add_argument('--permission', required=True, choices=['view', 'edit', 'owner'], 
                             help='Permission level')

    args = parser.parse_args()

    # Create application context
    with app.app_context():
        if args.command == 'assign':
            success = assign_permissions(args.email, args.project, args.permission)
        else:
            parser.print_help()
            sys.exit(1)
    
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main() 