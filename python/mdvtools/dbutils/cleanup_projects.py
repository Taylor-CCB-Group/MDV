#!/usr/bin/env python3
"""
Project Cleanup Script

This script provides utilities for cleaning up projects in the MDV database,
including removing deleted projects and handling projects with empty datasources.

As of writing, this is mostly intended for use in local development environments
- be very cautious about running in any kind of actual server deployment!
"""

import os
import sys
import argparse
import shutil
import json
from typing import List, Dict, Tuple
from datetime import datetime

from mdvtools.dbutils.mdv_server_app import app, is_valid_mdv_project
from mdvtools.dbutils.dbmodels import db, Project, File, UserProject
from mdvtools.logging_config import get_logger

logger = get_logger(__name__)


def list_deleted_projects(app) -> List[Project]:
    """
    Query and return all projects marked as deleted.
    
    Args:
        app: Flask application instance
        
    Returns:
        List of Project objects where is_deleted=True
    """
    with app.app_context():
        deleted_projects = Project.query.filter(Project.is_deleted == True).all()
        return deleted_projects


def list_projects_by_name_pattern(app, pattern: str) -> List[Project]:
    """
    Find all projects whose name matches the given pattern.
    
    Supports SQL LIKE wildcards (% for multiple chars, _ for single char).
    Also converts shell-style wildcards (* and ?) to SQL LIKE format.
    
    Args:
        app: Flask application instance
        pattern: Name pattern to match (e.g., "pbmc3k%", "*test*")
        
    Returns:
        List of Project objects matching the pattern
    """
    with app.app_context():
        # Convert shell-style wildcards to SQL LIKE wildcards
        sql_pattern = pattern.replace('*', '%').replace('?', '_')
        
        # Query projects matching the pattern
        projects = Project.query.filter(Project.name.like(sql_pattern)).all()
        logger.info(f"Found {len(projects)} projects matching pattern '{pattern}'")
        return projects


def list_orphaned_projects(app) -> List[Tuple[str, str]]:
    """
    Find all projects that exist in the filesystem but are not registered in the database.
    
    Args:
        app: Flask application instance
        
    Returns:
        List of tuples (project_path, project_name) for orphaned projects
    """
    orphaned = []
    
    with app.app_context():
        # Get the base directory for projects
        base_dir = app.config.get('projects_base_dir')
        if not base_dir or not os.path.exists(base_dir):
            logger.warning(f"Base directory not found or not configured: {base_dir}")
            return orphaned
        
        # Get all project paths from database
        db_projects = Project.query.all()
        db_paths = {project.path for project in db_projects}
        
        # Scan filesystem for project directories
        try:
            for item in os.listdir(base_dir):
                item_path = os.path.join(base_dir, item)
                
                # Check if it's a valid MDV project
                if is_valid_mdv_project(item_path):
                    # Check if this path is in the database
                    if item_path not in db_paths:
                        orphaned.append((item_path, item))
                        logger.info(f"Found orphaned project: {item_path}")
        except Exception as e:
            logger.error(f"Error scanning base directory {base_dir}: {e}")
    
    return orphaned


def list_test_generated_projects(app) -> List[Project]:
    """
    Find all projects that have test generation provenance in state.json.
    
    Args:
        app: Flask application instance
        
    Returns:
        List of Project objects with test generation provenance
    """
    test_projects = []
    
    with app.app_context():
        all_projects = Project.query.all()
        
        for project in all_projects:
            state_file = os.path.join(project.path, "state.json")
            
            if os.path.exists(state_file):
                try:
                    with open(state_file, 'r') as f:
                        state = json.load(f)
                        # Check if provenance metadata exists
                        if "provenance" in state:
                            test_projects.append(project)
                            logger.debug(f"Found test-generated project: {project.id} ({project.name})")
                except (json.JSONDecodeError, IOError) as e:
                    logger.warning(f"Error reading state.json for project {project.id} ({project.name}): {e}")
            else:
                logger.debug(f"No state.json found for project {project.id} at {project.path}")
    
    logger.info(f"Found {len(test_projects)} test-generated projects")
    return test_projects


def list_empty_datasource_projects(app) -> List[Project]:
    """
    Find all projects that have empty datasources (datasources.json contains []).
    
    Args:
        app: Flask application instance
        
    Returns:
        List of Project objects with empty datasources
    """
    with app.app_context():
        all_projects = Project.query.filter(Project.is_deleted == False).all()
        empty_projects = []
        
        for project in all_projects:
            datasources_file = os.path.join(project.path, "datasources.json")
            
            if os.path.exists(datasources_file):
                try:
                    with open(datasources_file, 'r') as f:
                        datasources = json.load(f)
                        # Check if datasources is an empty list
                        if isinstance(datasources, list) and len(datasources) == 0:
                            empty_projects.append(project)
                except (json.JSONDecodeError, IOError) as e:
                    logger.warning(f"Error reading datasources for project {project.id} ({project.name}): {e}")
            else:
                logger.debug(f"No datasources.json found for project {project.id} at {project.path}")
        
        return empty_projects


def delete_project_from_database(project_id: int, project_name: str) -> Tuple[bool, str]:
    """
    Delete a project and all related records from the database.
    Handles foreign key constraints by deleting related records first.
    
    Database Relationships:
    - File.project_id -> projects.id (NOT NULL foreign key)
    - UserProject.project_id -> projects.id (NOT NULL foreign key)
    
    These related records must be deleted BEFORE the project itself to avoid
    foreign key constraint violations.
    
    Args:
        project_id: Project ID to delete
        project_name: Project name (for logging)
        
    Returns:
        Tuple of (success: bool, error_message: str)
    """
    try:
        # Delete related File records first
        files = File.query.filter_by(project_id=project_id).all()
        if files:
            logger.info(f"  Deleting {len(files)} related file record(s)")
            for file in files:
                db.session.delete(file)
        
        # Delete related UserProject records
        user_projects = UserProject.query.filter_by(project_id=project_id).all()
        if user_projects:
            logger.info(f"  Deleting {len(user_projects)} related user-project record(s)")
            for user_project in user_projects:
                db.session.delete(user_project)
        
        # Now delete the project itself
        project = Project.query.get(project_id)
        if project:
            db.session.delete(project)
            db.session.commit()
            logger.info(f"  ✓ Removed from database (including {len(files)} files, {len(user_projects)} user associations)")
            return True, ""
        else:
            logger.warning(f"  ⚠ Project {project_id} not found in database")
            return False, "Project not found in database"
            
    except Exception as e:
        error_msg = f"Error removing from database: {e}"
        logger.error(f"  ✗ {error_msg}")
        db.session.rollback()
        return False, error_msg


def purge_deleted_projects(app, dry_run=False) -> Dict:
    """
    Remove deleted projects from both filesystem and database.
    
    Args:
        app: Flask application instance
        dry_run: If True, show what would be done without executing
        
    Returns:
        Dictionary with operation results
    """
    results = {
        'total': 0,
        'filesystem_removed': 0,
        'database_removed': 0,
        'errors': []
    }
    
    with app.app_context():
        deleted_projects = Project.query.filter(Project.is_deleted == True).all()
        results['total'] = len(deleted_projects)
        
        if dry_run:
            logger.info(f"[DRY RUN] Would purge {len(deleted_projects)} deleted projects")
            return results
        
        for project in deleted_projects:
            project_id = project.id
            project_name = project.name
            project_path = project.path
            
            logger.info(f"Processing project {project_id}: {project_name}")
            
            # Remove from filesystem
            if os.path.exists(project_path):
                try:
                    if os.path.isdir(project_path):
                        shutil.rmtree(project_path)
                        logger.info(f"  ✓ Removed directory: {project_path}")
                        results['filesystem_removed'] += 1
                    else:
                        os.remove(project_path)
                        logger.info(f"  ✓ Removed file: {project_path}")
                        results['filesystem_removed'] += 1
                except Exception as e:
                    error_msg = f"Error removing filesystem path: {e}"
                    logger.error(f"  ✗ {error_msg}")
                    results['errors'].append((project_id, project_name, error_msg))
            else:
                logger.info(f"  ⚠ Path does not exist: {project_path}")
            
            # Remove from database (including related records)
            success, error_msg = delete_project_from_database(project_id, project_name)
            if success:
                results['database_removed'] += 1
            else:
                results['errors'].append((project_id, project_name, error_msg))
    
    return results


def purge_orphaned_projects(app, dry_run=False) -> Dict:
    """
    Remove orphaned project directories (exist in filesystem but not in database).
    
    Args:
        app: Flask application instance
        dry_run: If True, show what would be done without executing
        
    Returns:
        Dictionary with operation results
    """
    results = {
        'total': 0,
        'removed': 0,
        'errors': []
    }
    
    with app.app_context():
        orphaned_projects = list_orphaned_projects(app)
        results['total'] = len(orphaned_projects)
        
        if dry_run:
            logger.info(f"[DRY RUN] Would remove {len(orphaned_projects)} orphaned project directories")
            return results
        
        for project_path, project_name in orphaned_projects:
            logger.info(f"Processing orphaned project: {project_name} at {project_path}")
            
            if os.path.exists(project_path):
                try:
                    shutil.rmtree(project_path)
                    logger.info(f"  ✓ Removed directory: {project_path}")
                    results['removed'] += 1
                except Exception as e:
                    error_msg = f"Error removing orphaned project: {e}"
                    logger.error(f"  ✗ {error_msg}")
                    results['errors'].append((project_name, project_name, error_msg))
            else:
                logger.warning(f"  ⚠ Path no longer exists: {project_path}")
    
    return results


def purge_empty_datasource_projects(app, hard_delete=False, dry_run=False) -> Dict:
    """
    Handle projects with empty datasources.
    
    Args:
        app: Flask application instance
        hard_delete: If True, permanently delete. If False, soft-delete.
        dry_run: If True, show what would be done without executing
        
    Returns:
        Dictionary with operation results
    """
    results = {
        'total': 0,
        'soft_deleted': 0,
        'hard_deleted': 0,
        'errors': []
    }
    
    with app.app_context():
        empty_projects = list_empty_datasource_projects(app)
        results['total'] = len(empty_projects)
        
        action = "hard-delete" if hard_delete else "soft-delete"
        if dry_run:
            logger.info(f"[DRY RUN] Would {action} {len(empty_projects)} projects with empty datasources")
            return results
        
        for project in empty_projects:
            project_id = project.id
            project_name = project.name
            project_path = project.path
            
            logger.info(f"Processing project {project_id}: {project_name}")
            
            if hard_delete:
                # Remove from filesystem
                if os.path.exists(project_path):
                    try:
                        if os.path.isdir(project_path):
                            shutil.rmtree(project_path)
                            logger.info(f"  ✓ Removed directory: {project_path}")
                        else:
                            os.remove(project_path)
                            logger.info(f"  ✓ Removed file: {project_path}")
                    except Exception as e:
                        error_msg = f"Error removing filesystem path: {e}"
                        logger.error(f"  ✗ {error_msg}")
                        results['errors'].append((project_id, project_name, error_msg))
                        continue
                
                # Remove from database (including related records)
                success, error_msg = delete_project_from_database(project_id, project_name)
                if success:
                    results['hard_deleted'] += 1
                else:
                    results['errors'].append((project_id, project_name, error_msg))
            else:
                # Soft delete
                try:
                    project.is_deleted = True
                    project.deleted_timestamp = datetime.now()
                    db.session.commit()
                    logger.info(f"  ✓ Soft-deleted (marked as deleted in database)")
                    results['soft_deleted'] += 1
                except Exception as e:
                    error_msg = f"Error soft-deleting project {project_id} ({project_name}): {e}"
                    logger.error(f"  ✗ {error_msg}")
                    results['errors'].append((project_id, project_name, error_msg))
                    db.session.rollback()
    
    return results


def print_project_list(projects: List[Project], title: str):
    """
    Print a formatted list of projects.
    
    Args:
        projects: List of Project objects to display
        title: Title for the list
    """
    print(f"\n{'='*80}")
    print(f"{title}")
    print(f"{'='*80}")
    print(f"Found {len(projects)} project(s)\n")
    
    if not projects:
        print("No projects found.")
        return
    
    for i, project in enumerate(projects, 1):
        exists = "✓ EXISTS" if os.path.exists(project.path) else "✗ NOT FOUND"
        print(f"{i}. [{exists}] ID {project.id}: {project.name}")
        print(f"   Path: {project.path}")
        
        if project.is_deleted and project.deleted_timestamp:
            print(f"   Deleted: {project.deleted_timestamp}")
        
        # Try to get datasource count
        datasources_file = os.path.join(project.path, "datasources.json")
        if os.path.exists(datasources_file):
            try:
                with open(datasources_file, 'r') as f:
                    datasources = json.load(f)
                    if isinstance(datasources, list):
                        print(f"   Datasources: {len(datasources)}")
            except:
                print(f"   Datasources: Error reading file")
        
        # Try to get provenance metadata from state.json
        state_file = os.path.join(project.path, "state.json")
        if os.path.exists(state_file):
            try:
                with open(state_file, 'r') as f:
                    state = json.load(f)
                    if "provenance" in state:
                        prov = state["provenance"]
                        # Format provenance info compactly
                        source = prov.get("source", "unknown")
                        created_by = prov.get("created_by", "unknown")
                        created_at = prov.get("created_at", "unknown")
                        
                        # Build parameter string
                        params = prov.get("parameters", {})
                        param_parts = []
                        if "dataset" in params and params["dataset"]:
                            param_parts.append(f"dataset={params['dataset']}")
                        if "n_cells" in params:
                            param_parts.append(f"{params['n_cells']:,} cells")
                        if "n_genes" in params:
                            param_parts.append(f"{params['n_genes']:,} genes")
                        
                        param_str = ", ".join(param_parts) if param_parts else ""
                        prov_str = f"{created_by} ({source}"
                        if param_str:
                            prov_str += f", {param_str}"
                        prov_str += ")"
                        
                        print(f"   Provenance: {prov_str}")
                        if created_at != "unknown":
                            # Try to parse and format datetime
                            try:
                                from datetime import datetime as dt
                                created_dt = dt.fromisoformat(created_at.replace('Z', '+00:00'))
                                print(f"   Created: {created_dt.strftime('%Y-%m-%d %H:%M:%S')}")
                            except:
                                print(f"   Created: {created_at}")
            except (json.JSONDecodeError, IOError):
                pass  # Silently skip if can't read state.json
        
        print()


def print_orphaned_project_list(orphaned_projects: List[Tuple[str, str]], title: str):
    """
    Print a formatted list of orphaned projects.
    
    Args:
        orphaned_projects: List of (path, name) tuples
        title: Title for the list
    """
    print(f"\n{'='*80}")
    print(f"{title}")
    print(f"{'='*80}")
    print(f"Found {len(orphaned_projects)} orphaned project(s)\n")
    
    if not orphaned_projects:
        print("No orphaned projects found.")
        return
    
    for i, (project_path, project_name) in enumerate(orphaned_projects, 1):
        print(f"{i}. {project_name}")
        print(f"   Path: {project_path}")
        
        # Try to get datasource count
        datasources_file = os.path.join(project_path, "datasources.json")
        if os.path.exists(datasources_file):
            try:
                with open(datasources_file, 'r') as f:
                    datasources = json.load(f)
                    if isinstance(datasources, list):
                        print(f"   Datasources: {len(datasources)}")
            except:
                print(f"   Datasources: Error reading file")
        
        print()


def print_results(results: Dict, operation: str):
    """
    Print operation results summary.
    
    Args:
        results: Dictionary with operation results
        operation: Name of the operation
    """
    print(f"\n{'='*80}")
    print(f"SUMMARY: {operation}")
    print(f"{'='*80}")
    print(f"Total projects processed: {results['total']}")
    
    if 'filesystem_removed' in results:
        print(f"Filesystem removals: {results['filesystem_removed']}")
    if 'database_removed' in results:
        print(f"Database removals: {results['database_removed']}")
    if 'removed' in results:
        print(f"Removed: {results['removed']}")
    if 'soft_deleted' in results:
        print(f"Soft-deleted: {results['soft_deleted']}")
    if 'hard_deleted' in results:
        print(f"Hard-deleted: {results['hard_deleted']}")
    
    if results['errors']:
        print(f"\nErrors encountered: {len(results['errors'])}")
        # Handle both 2-tuple and 3-tuple error formats
        for error_info in results['errors']:
            if len(error_info) == 3:
                project_id, project_name, error = error_info
                print(f"  - Project {project_id} ({project_name}): {error}")
            elif len(error_info) == 2:
                project_name, error = error_info
                print(f"  - {project_name}: {error}")
            else:
                print(f"  - {error_info}")
    else:
        print("\n✓ All operations completed successfully")


def confirm_action(message: str) -> bool:
    """
    Prompt user for confirmation.
    
    Args:
        message: Confirmation message to display
        
    Returns:
        True if user confirms, False otherwise
    """
    response = input(f"\n{message} (yes/no): ").strip().lower()
    return response in ['yes', 'y']


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description='Cleanup utility for MDV projects database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List deleted projects
  %(prog)s --list-deleted
  
  # List orphaned projects (in filesystem but not database)
  %(prog)s --list-orphaned
  
  # Preview what would be purged
  %(prog)s --purge-deleted --dry-run
  
  # Actually purge deleted projects
  %(prog)s --purge-deleted --confirm
  
  # List projects with empty datasources
  %(prog)s --list-empty
  
  # Soft-delete projects with empty datasources
  %(prog)s --purge-empty-datasources
  
  # Hard-delete (purge) projects with empty datasources
  %(prog)s --purge-empty-datasources --hard-delete --confirm
  
  # Remove orphaned project directories
  %(prog)s --purge-orphaned --confirm
        """
    )
    
    # Operation modes
    parser.add_argument('--list-deleted', action='store_true',
                        help='List all projects marked as deleted')
    parser.add_argument('--list-empty', action='store_true',
                        help='List projects with empty datasources')
    parser.add_argument('--list-orphaned', action='store_true',
                        help='List projects that exist in filesystem but not in database')
    parser.add_argument('--list-test-generated', action='store_true',
                        help='List projects with test generation provenance metadata')
    parser.add_argument('--purge-deleted', action='store_true',
                        help='Remove deleted projects from filesystem and database')
    parser.add_argument('--purge-empty-datasources', action='store_true',
                        help='Delete projects with empty datasources')
    parser.add_argument('--purge-orphaned', action='store_true',
                        help='Remove orphaned project directories from filesystem')
    parser.add_argument('--delete-and-purge', action='store_true',
                        help='Delete and purge projects in one step (use with --filter-name or --list-test-generated)')
    
    # Options
    parser.add_argument('--hard-delete', action='store_true',
                        help='Use with --purge-empty-datasources for immediate purge (default is soft-delete)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be done without executing')
    parser.add_argument('--confirm', action='store_true',
                        help='Skip confirmation prompts (use with caution)')
    
    # Filters
    parser.add_argument('--filter-name', type=str, metavar='PATTERN',
                        help='Filter projects by name pattern (supports * and %% wildcards)')
    
    args = parser.parse_args()
    
    # Ensure app is available
    if app is None:
        logger.error("Flask app is not initialized. Cannot proceed.")
        sys.exit(1)
    
    # Check that at least one operation is specified
    if not any([args.list_deleted, args.list_empty, args.list_orphaned, args.list_test_generated,
                args.purge_deleted, args.purge_empty_datasources, args.purge_orphaned, args.delete_and_purge]):
        parser.print_help()
        sys.exit(1)
    
    # Validate delete-and-purge requires a filter
    if args.delete_and_purge and not (args.filter_name or args.list_test_generated):
        print("Error: --delete-and-purge requires either --filter-name or --list-test-generated for safety")
        print("This prevents accidentally deleting all projects.")
        sys.exit(1)
    
    # Execute operations
    try:
        # Helper function to apply name filter
        def apply_name_filter(projects: List[Project]) -> List[Project]:
            if args.filter_name:
                name_filtered = list_projects_by_name_pattern(app, args.filter_name)
                # Get intersection of projects
                name_filtered_ids = {p.id for p in name_filtered}
                return [p for p in projects if p.id in name_filtered_ids]
            return projects
        
        if args.list_deleted:
            projects = list_deleted_projects(app)
            projects = apply_name_filter(projects)
            title = "DELETED PROJECTS" + (f" (filtered by name: {args.filter_name})" if args.filter_name else "")
            print_project_list(projects, title)
        
        if args.list_empty:
            projects = list_empty_datasource_projects(app)
            projects = apply_name_filter(projects)
            title = "PROJECTS WITH EMPTY DATASOURCES" + (f" (filtered by name: {args.filter_name})" if args.filter_name else "")
            print_project_list(projects, title)
        
        if args.list_orphaned:
            orphaned = list_orphaned_projects(app)
            # For orphaned projects, we can't filter by name in the same way since they're tuples
            # but we can filter by the name part of the tuple
            if args.filter_name:
                sql_pattern = args.filter_name.replace('*', '%').replace('?', '_')
                # Convert SQL LIKE pattern to Python regex for simple matching
                import re
                regex_pattern = sql_pattern.replace('%', '.*').replace('_', '.')
                orphaned = [(path, name) for path, name in orphaned if re.match(regex_pattern, name)]
                title = f"ORPHANED PROJECTS (Filesystem Only) (filtered by name: {args.filter_name})"
            else:
                title = "ORPHANED PROJECTS (Filesystem Only)"
            print_orphaned_project_list(orphaned, title)
        
        if args.list_test_generated:
            projects = list_test_generated_projects(app)
            projects = apply_name_filter(projects)
            title = "TEST-GENERATED PROJECTS" + (f" (filtered by name: {args.filter_name})" if args.filter_name else "")
            print_project_list(projects, title)
        
        if args.purge_deleted:
            # Apply name filter if specified
            if args.filter_name:
                logger.info(f"Applying name filter: {args.filter_name}")
            
            if not args.dry_run and not args.confirm:
                filter_msg = f" matching '{args.filter_name}'" if args.filter_name else ""
                if not confirm_action(f"⚠️  This will PERMANENTLY DELETE{filter_msg} projects from filesystem and database. Continue?"):
                    print("Operation cancelled.")
                    sys.exit(0)
            
            # Get deleted projects and apply filter
            with app.app_context():
                projects_to_purge = list_deleted_projects(app)
                projects_to_purge = apply_name_filter(projects_to_purge)
                
                # Create a modified version of purge that works with filtered list
                results = {
                    'total': len(projects_to_purge),
                    'filesystem_removed': 0,
                    'database_removed': 0,
                    'errors': []
                }
                
                if args.dry_run:
                    logger.info(f"[DRY RUN] Would purge {len(projects_to_purge)} deleted projects")
                else:
                    for project in projects_to_purge:
                        project_id = project.id
                        project_name = project.name
                        project_path = project.path
                        
                        logger.info(f"Processing project {project_id}: {project_name}")
                        
                        # Remove from filesystem
                        if os.path.exists(project_path):
                            try:
                                if os.path.isdir(project_path):
                                    shutil.rmtree(project_path)
                                    logger.info(f"  ✓ Removed directory: {project_path}")
                                    results['filesystem_removed'] += 1
                                else:
                                    os.remove(project_path)
                                    logger.info(f"  ✓ Removed file: {project_path}")
                                    results['filesystem_removed'] += 1
                            except Exception as e:
                                error_msg = f"Error removing filesystem path: {e}"
                                logger.error(f"  ✗ {error_msg}")
                                results['errors'].append((project_id, project_name, error_msg))
                        else:
                            logger.info(f"  ⚠ Path does not exist: {project_path}")
                        
                        # Remove from database (including related records)
                        success, error_msg = delete_project_from_database(project_id, project_name)
                        if success:
                            results['database_removed'] += 1
                        else:
                            results['errors'].append((project_id, project_name, error_msg))
            
            filter_title = f" (filtered by name: {args.filter_name})" if args.filter_name else ""
            print_results(results, f"PURGE DELETED PROJECTS{filter_title}")
        
        if args.purge_empty_datasources:
            action = "HARD-DELETE" if args.hard_delete else "SOFT-DELETE"
            if not args.dry_run and not args.confirm:
                message = f"⚠️  This will {action} projects with empty datasources. Continue?"
                if not confirm_action(message):
                    print("Operation cancelled.")
                    sys.exit(0)
            
            results = purge_empty_datasource_projects(app, hard_delete=args.hard_delete, dry_run=args.dry_run)
            print_results(results, f"{action} EMPTY DATASOURCE PROJECTS")
        
        if args.purge_orphaned:
            if not args.dry_run and not args.confirm:
                if not confirm_action("⚠️  This will PERMANENTLY DELETE orphaned project directories. Continue?"):
                    print("Operation cancelled.")
                    sys.exit(0)
            
            results = purge_orphaned_projects(app, dry_run=args.dry_run)
            print_results(results, "PURGE ORPHANED PROJECTS")
        
        if args.delete_and_purge:
            # Get projects matching the criteria
            with app.app_context():
                if args.list_test_generated:
                    projects_to_delete = list_test_generated_projects(app)
                    criteria_msg = "test-generated projects"
                else:
                    # Must have filter_name due to validation above
                    projects_to_delete = list_projects_by_name_pattern(app, args.filter_name)
                    criteria_msg = f"projects matching '{args.filter_name}'"
                
                # Apply additional name filter if both are specified
                if args.filter_name and args.list_test_generated:
                    projects_to_delete = apply_name_filter(projects_to_delete)
                    criteria_msg = f"test-generated projects matching '{args.filter_name}'"
                
                # Confirm action
                if not args.dry_run and not args.confirm:
                    msg = f"⚠️  This will DELETE and PURGE {len(projects_to_delete)} {criteria_msg}. Continue?"
                    if not confirm_action(msg):
                        print("Operation cancelled.")
                        sys.exit(0)
                
                # Execute delete and purge
                results = {
                    'total': len(projects_to_delete),
                    'filesystem_removed': 0,
                    'database_removed': 0,
                    'errors': []
                }
                
                if args.dry_run:
                    logger.info(f"[DRY RUN] Would delete and purge {len(projects_to_delete)} {criteria_msg}")
                else:
                    for project in projects_to_delete:
                        project_id = project.id
                        project_name = project.name
                        project_path = project.path
                        
                        logger.info(f"Processing project {project_id}: {project_name}")
                        
                        # Remove from filesystem
                        if os.path.exists(project_path):
                            try:
                                if os.path.isdir(project_path):
                                    shutil.rmtree(project_path)
                                    logger.info(f"  ✓ Removed directory: {project_path}")
                                    results['filesystem_removed'] += 1
                                else:
                                    os.remove(project_path)
                                    logger.info(f"  ✓ Removed file: {project_path}")
                                    results['filesystem_removed'] += 1
                            except Exception as e:
                                error_msg = f"Error removing filesystem path: {e}"
                                logger.error(f"  ✗ {error_msg}")
                                results['errors'].append((project_id, project_name, error_msg))
                        else:
                            logger.info(f"  ⚠ Path does not exist: {project_path}")
                        
                        # Remove from database (including related records)
                        success, error_msg = delete_project_from_database(project_id, project_name)
                        if success:
                            results['database_removed'] += 1
                        else:
                            results['errors'].append((project_id, project_name, error_msg))
                
                print_results(results, f"DELETE AND PURGE: {criteria_msg}")
    
    except Exception as e:
        logger.exception(f"Fatal error during cleanup operation: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()

