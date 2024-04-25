import os
import json
from mdvtools.mdvproject import MDVProject
from mdvtools.dbutils.dbmodels import Project

def load_projects_from_config(base_dir):
    config_file = os.path.join(base_dir, 'project_config.json')
    config_projects = []

    try:
        if os.path.exists(config_file):
            with open(config_file, 'r') as f:
                config = json.load(f)
                for d in config.get('projects', []):
                    print(f"Adding project '{d}' from config file")
                    try:
                        config_projects.append(MDVProject(os.path.join(base_dir, d)))
                    except Exception as e:
                        print(f"Error adding project '{d}' from config file: {e}")
        else:
            print("Config file doesn't exist.")
            with open(config_file, 'w') as f:
                json.dump({'projects': []}, f)
    except Exception as e:
        print(f"Error loading projects from config: {e}")

    return config_projects

def create_projects_from_db(base_dir):
    db_projects = []

    try:
        # Assuming Project.query.all() returns a list of project names
        project_names = [project.name for project in Project.query.all()]

        for project_name in project_names:
            project_path = os.path.join(base_dir, project_name)
            try:
                if os.path.exists(project_path):
                    db_projects.append(MDVProject(project_path))
                else:
                    print(f"Error: Project path '{project_path}' does not exist.")
            except Exception as e:
                print(f"Error creating project '{project_name}': {e}")
    except Exception as e:
        print(f"Error retrieving projects from database: {e}")

    return db_projects


def create_all_projects(base_dir):
    try:
        db_projects = create_projects_from_db(base_dir)
        config_projects = load_projects_from_config(base_dir)
        return db_projects + config_projects
    except Exception as e:
        print(f"Error creating all projects: {e}")
        return []

