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
from typing import List, Dict, Tuple, Optional, Any, Union
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


def delete_project_from_database(app, project_id: int, project_name: str, dry_run: bool = False) -> Dict[str, Any]:
    """
    Delete a project and all related records from the database.
    Handles foreign key constraints by deleting related records first.
    
    Database Relationships:
    - File.project_id -> projects.id (NOT NULL foreign key)
    - UserProject.project_id -> projects.id (NOT NULL foreign key)
    
    These related records must be deleted BEFORE the project itself to avoid
    foreign key constraint violations.
    
    Args:
        app: Flask application instance
        project_id: Project ID to delete
        project_name: Project name (for logging)
        dry_run: If True, only preview without executing
        
    Returns:
        Dictionary with success count, failed count, and errors
    """
    results = {
        'success_count': 0,
        'failed_count': 0,
        'errors': []
    }
    
    try:
        with app.app_context():
            if dry_run:
                # Just log what would be done
                files = File.query.filter_by(project_id=project_id).all()
                user_projects = UserProject.query.filter_by(project_id=project_id).all()
                logger.info(f"  [DRY RUN] Would delete {len(files)} file record(s), {len(user_projects)} user-project record(s), and project {project_id}")
                results['success_count'] = 1
            else:
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
                    results['success_count'] = 1
                else:
                    logger.warning(f"  ⚠ Project {project_id} not found in database")
                    results['failed_count'] = 1
                    results['errors'].append((project_name, "Project not found in database"))
            
    except Exception as e:
        error_msg = f"Error removing from database: {e}"
        logger.error(f"  ✗ {error_msg}")
        results['failed_count'] = 1
        results['errors'].append((project_name, str(e)))
        if not dry_run:
            db.session.rollback()
    
    return results


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
            delete_result = delete_project_from_database(app, project_id, project_name, dry_run=False)
            if delete_result['success_count'] > 0:
                results['database_removed'] += 1
            else:
                if delete_result['errors']:
                    results['errors'].extend(delete_result['errors'])
    
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
                delete_result = delete_project_from_database(app, project_id, project_name, dry_run=False)
                if delete_result['success_count'] > 0:
                    results['hard_deleted'] += 1
                else:
                    if delete_result['errors']:
                        results['errors'].extend(delete_result['errors'])
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
    
    # Handle both old and new result formats
    if 'total' in results:
        print(f"Total projects processed: {results['total']}")
    else:
        total = results.get('success_count', 0) + results.get('failed_count', 0)
        print(f"Total projects processed: {total}")
    
    # New format
    if 'success_count' in results:
        print(f"Successful: {results['success_count']}")
    if 'failed_count' in results and results['failed_count'] > 0:
        print(f"Failed: {results['failed_count']}")
    
    # Old format (for backwards compatibility)
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
    
    if results.get('errors'):
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


def apply_name_filter(projects: List[Project], pattern: str) -> List[Project]:
    """
    Filter projects by name pattern.
    
    Args:
        projects: List of Project objects
        pattern: Name pattern (supports * and % wildcards)
        
    Returns:
        Filtered list of projects
    """
    # Convert wildcards to SQL LIKE pattern
    sql_pattern = pattern.replace('*', '%')
    
    # Filter projects
    filtered = []
    for project in projects:
        # Use SQL LIKE pattern matching
        if '%' in sql_pattern:
            # Simple pattern matching
            import re
            regex_pattern = sql_pattern.replace('%', '.*')
            if re.match(regex_pattern, project.name):
                filtered.append(project)
        else:
            # Exact match
            if project.name == sql_pattern:
                filtered.append(project)
    
    return filtered


def select_projects(app, selector: str, name_filter: Optional[str] = None) -> Union[List[Project], Tuple[List[Project], List[Tuple[str, str]]]]:
    """
    Universal project selector.
    
    Args:
        app: Flask application instance
        selector: Selector type (deleted, empty, orphaned, test-generated, all)
        name_filter: Optional name pattern filter
        
    Returns:
        For orphaned selector: tuple (empty list, orphaned_list)
        For all other selectors: List of Project objects
    """
    logger.info(f"Selecting projects with selector: {selector}, name_filter: {name_filter}")
    
    if selector == 'deleted':
        projects = list_deleted_projects(app)
    elif selector == 'empty':
        projects = list_empty_datasource_projects(app)
    elif selector == 'orphaned':
        # Orphaned projects are special - they return tuples (path, name)
        orphaned = list_orphaned_projects(app)
        if name_filter:
            # Filter orphaned projects by name
            sql_pattern = name_filter.replace('*', '%')
            import re
            regex_pattern = sql_pattern.replace('%', '.*')
            filtered_orphaned = []
            for path, name in orphaned:
                if '%' in sql_pattern:
                    if re.match(regex_pattern, name):
                        filtered_orphaned.append((path, name))
                else:
                    if name == sql_pattern:
                        filtered_orphaned.append((path, name))
            return ([], filtered_orphaned)  # Return empty projects list and filtered orphaned
        return ([], orphaned)  # Return empty projects list and all orphaned
    elif selector == 'test-generated':
        projects = list_test_generated_projects(app)
    elif selector == 'all':
        # Get all projects
        with app.app_context():
            projects = Project.query.all()
    else:
        logger.error(f"Unknown selector: {selector}")
        return []
    
    # Apply name filter if provided (not for orphaned, handled above)
    if name_filter and selector != 'orphaned':
        projects = apply_name_filter(projects, name_filter)
    
    return projects


def action_list(app, projects: List[Project], orphaned_projects: List[tuple], selector: str):
    """
    List projects.
    
    Args:
        app: Flask application instance
        projects: List of Project objects
        orphaned_projects: List of orphaned project tuples
        selector: Selector type for title
    """
    if selector == 'orphaned':
        title = f"{selector.upper()} PROJECTS"
        print_orphaned_project_list(orphaned_projects, title)
    else:
        title = f"{selector.upper()} PROJECTS"
        print_project_list(projects, title)


def action_soft_delete(app, projects: List[Project], dry_run: bool = False) -> Dict[str, Any]:
    """
    Soft-delete projects (mark as deleted).
    
    Args:
        app: Flask application instance
        projects: List of Project objects
        dry_run: If True, only preview without executing
        
    Returns:
        Dictionary with success count, failed count, and errors
    """
    from datetime import datetime
    
    results = {
        'success_count': 0,
        'failed_count': 0,
        'errors': []
    }
    
    with app.app_context():
        for project in projects:
            try:
                if dry_run:
                    print(f"Would mark as deleted: {project.name} (ID: {project.id}, Path: {project.path})")
                    results['success_count'] += 1
                else:
                    project.is_deleted = True
                    project.deleted_timestamp = datetime.now()
                    db.session.commit()
                    logger.info(f"Marked project as deleted: {project.name} (ID: {project.id})")
                    print(f"✓ Marked as deleted: {project.name}")
                    results['success_count'] += 1
            except Exception as e:
                error_msg = f"Failed to mark {project.name} as deleted: {str(e)}"
                logger.error(error_msg)
                results['failed_count'] += 1
                results['errors'].append((project.name, str(e)))
                if not dry_run:
                    db.session.rollback()
    
    return results


def action_purge(app, projects: List[Project], orphaned_projects: List[tuple], 
                 selector: str, dry_run: bool = False) -> Dict[str, Any]:
    """
    Purge projects completely (filesystem + database).
    
    Args:
        app: Flask application instance
        projects: List of Project objects
        orphaned_projects: List of orphaned project tuples
        selector: Selector type to determine purge strategy
        dry_run: If True, only preview without executing
        
    Returns:
        Dictionary with success count, failed count, and errors
    """
    results = {
        'success_count': 0,
        'failed_count': 0,
        'errors': []
    }
    
    # Handle orphaned projects separately
    if selector == 'orphaned':
        with app.app_context():
            for project_path, project_name in orphaned_projects:
                try:
                    if dry_run:
                        print(f"Would remove orphaned directory: {project_name} ({project_path})")
                        results['success_count'] += 1
                    else:
                        if os.path.isdir(project_path):
                            shutil.rmtree(project_path)
                            logger.info(f"Removed orphaned directory: {project_path}")
                            print(f"✓ Removed orphaned directory: {project_name}")
                            results['success_count'] += 1
                        else:
                            logger.warning(f"Orphaned path no longer exists: {project_path}")
                            results['success_count'] += 1
                except Exception as e:
                    error_msg = f"Failed to remove orphaned directory {project_name}: {str(e)}"
                    logger.error(error_msg)
                    results['failed_count'] += 1
                    results['errors'].append((project_name, str(e)))
        return results
    
    # Handle database projects
    with app.app_context():
        for project in projects:
            try:
                if dry_run:
                    print(f"Would purge: {project.name} (ID: {project.id}, Path: {project.path})")
                    results['success_count'] += 1
                else:
                    # Remove from filesystem
                    if os.path.exists(project.path):
                        if os.path.isdir(project.path):
                            shutil.rmtree(project.path)
                        else:
                            os.remove(project.path)
                        logger.info(f"Removed from filesystem: {project.path}")
                    
                    # Delete from database (handles foreign keys)
                    delete_result = delete_project_from_database(app, project.id, project.name, dry_run=False)
                    
                    if delete_result['success_count'] > 0:
                        print(f"✓ Purged: {project.name}")
                        results['success_count'] += 1
                    else:
                        results['failed_count'] += 1
                        if delete_result['errors']:
                            results['errors'].extend(delete_result['errors'])
                        
            except Exception as e:
                error_msg = f"Failed to purge {project.name}: {str(e)}"
                logger.error(error_msg)
                results['failed_count'] += 1
                results['errors'].append((project.name, str(e)))
                if not dry_run:
                    db.session.rollback()
    
    return results


def action_restore(app, projects: List[Project], dry_run: bool = False) -> Dict[str, Any]:
    """
    Restore deleted projects (undelete).
    
    Args:
        app: Flask application instance
        projects: List of Project objects
        dry_run: If True, only preview without executing
        
    Returns:
        Dictionary with success count, failed count, and errors
    """
    results = {
        'success_count': 0,
        'failed_count': 0,
        'errors': []
    }
    
    with app.app_context():
        for project in projects:
            try:
                if dry_run:
                    print(f"Would restore: {project.name} (ID: {project.id}, Path: {project.path})")
                    results['success_count'] += 1
                else:
                    project.is_deleted = False
                    project.deleted_timestamp = None
                    db.session.commit()
                    logger.info(f"Restored project: {project.name} (ID: {project.id})")
                    print(f"✓ Restored: {project.name}")
                    results['success_count'] += 1
            except Exception as e:
                error_msg = f"Failed to restore {project.name}: {str(e)}"
                logger.error(error_msg)
                results['failed_count'] += 1
                results['errors'].append((project.name, str(e)))
                if not dry_run:
                    db.session.rollback()
    
    return results


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description='Cleanup utility for MDV projects database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List projects
  %(prog)s --list deleted
  %(prog)s --list empty --filter-name "pbmc*"
  %(prog)s --list test-generated
  
  # Soft-delete (mark as deleted)
  %(prog)s --soft-delete empty --dry-run
  %(prog)s --soft-delete test-generated --confirm
  
  # Purge (complete removal)
  %(prog)s --purge deleted --confirm
  %(prog)s --purge empty --filter-name "pbmc*" --confirm
  
  # Restore (undelete)
  %(prog)s --restore deleted --filter-name "important*"
        """
    )
    
    # Actions (mutually exclusive)
    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument('--list', 
                             choices=['deleted', 'empty', 'orphaned', 'test-generated', 'all'],
                             metavar='SELECTOR',
                             help='List projects (choices: deleted, empty, orphaned, test-generated, all)')
    action_group.add_argument('--soft-delete',
                             choices=['empty', 'test-generated', 'all'],
                             metavar='SELECTOR',
                             help='Mark projects as deleted (choices: empty, test-generated, all)')
    action_group.add_argument('--purge',
                             choices=['deleted', 'empty', 'orphaned', 'test-generated', 'all'],
                             metavar='SELECTOR',
                             help='Remove projects completely (choices: deleted, empty, orphaned, test-generated, all)')
    action_group.add_argument('--restore',
                             choices=['deleted', 'all'],
                             metavar='SELECTOR',
                             help='Restore deleted projects (choices: deleted, all)')
    
    # Options
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
    
    # Determine action and selector
    if args.list:
        action = 'list'
        selector = args.list
    elif args.soft_delete:
        action = 'soft_delete'
        selector = args.soft_delete
    elif args.purge:
        action = 'purge'
        selector = args.purge
    elif args.restore:
        action = 'restore'
        selector = args.restore
    else:
        parser.print_help()
        sys.exit(1)
    
    # Safety check for 'all' selector
    if selector == 'all' and action != 'list':
        if not args.confirm and not args.dry_run:
            print(f"Error: --{action.replace('_', '-')} all requires --confirm or --dry-run flag for safety")
            print("This prevents accidentally affecting all projects.")
            sys.exit(1)
    
    # Execute operations
    try:
        # Select projects based on selector and filter
        result = select_projects(app, selector, args.filter_name)
        
        # Handle orphaned projects separately (they return a tuple)
        projects: List[Project]
        orphaned_projects: List[Tuple[str, str]]
        
        if selector == 'orphaned':
            # Type narrowing: result is Tuple[List[Project], List[Tuple[str, str]]]
            assert isinstance(result, tuple)
            projects, orphaned_projects = result
        else:
            # Type narrowing: result is List[Project]
            assert isinstance(result, list)
            projects = result
            orphaned_projects = []
        
        # Print selection summary
        if args.filter_name:
            filter_msg = f" (filtered by name: {args.filter_name})"
        else:
            filter_msg = ""
        
        if selector == 'orphaned':
            count_msg = f"Found {len(orphaned_projects)} {selector} project(s){filter_msg}"
        else:
            count_msg = f"Found {len(projects)} {selector} project(s){filter_msg}"
        
        logger.info(count_msg)
        print(f"\n{count_msg}")
        
        # Execute action
        if action == 'list':
            action_list(app, projects, orphaned_projects, selector)
        
        elif action == 'soft_delete':
            # Confirm if needed
            if not args.dry_run and not args.confirm:
                if not confirm_action(f"⚠️  Mark {len(projects)} project(s) as deleted?"):
                    print("Operation cancelled.")
                    sys.exit(0)
            
            results = action_soft_delete(app, projects, args.dry_run)
            print_results(results, f"SOFT-DELETE {selector.upper()} PROJECTS{filter_msg}")
        
        elif action == 'purge':
            # Confirm if needed
            if not args.dry_run and not args.confirm:
                total = len(orphaned_projects) if selector == 'orphaned' else len(projects)
                if not confirm_action(f"⚠️  PERMANENTLY DELETE {total} project(s) from filesystem and database?"):
                    print("Operation cancelled.")
                    sys.exit(0)
            
            results = action_purge(app, projects, orphaned_projects, selector, args.dry_run)
            print_results(results, f"PURGE {selector.upper()} PROJECTS{filter_msg}")
        
        elif action == 'restore':
            # Confirm if needed
            if not args.dry_run and not args.confirm:
                if not confirm_action(f"⚠️  Restore {len(projects)} project(s)?"):
                    print("Operation cancelled.")
                    sys.exit(0)
            
            results = action_restore(app, projects, args.dry_run)
            print_results(results, f"RESTORE {selector.upper()} PROJECTS{filter_msg}")
    
    except Exception as e:
        logger.exception(f"Fatal error during cleanup operation: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()

