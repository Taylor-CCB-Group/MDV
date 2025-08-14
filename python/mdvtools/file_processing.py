import os
import scanpy as sc
from mdvtools.conversions import convert_scanpy_to_mdv

def datasource_processing(project, filepath, original_filename, view, replace, supplied_only, resume=False):
    """
    Enhanced datasource processing with resumption capability.
    
    Args:
        project (Any): The project instance to which the datasource will be added.
        filepath (str): Path to the CSV file.
        original_filename (str): Original name of the uploaded file.
        view (str): The view to associate the datasource with.
        replace (bool): Whether to replace existing data.
        supplied_only (bool): Whether to use only supplied columns.
        resume (bool): Whether to attempt resuming an interrupted upload.

    Returns:
        dict: A dictionary containing 'success' and metadata on success.
        
    Raises:
        ValidationError: If any exception occurs during the processing, encapsulating
                         the error message with a 400 status code.
    """
    print(f"Processing datasource: {original_filename} in view: {view} (resume={resume})")
    print(f"Filepath: {filepath}")
    print(f"Replace: {replace}, Supplied only: {supplied_only}")

    try:
        # Check for existing progress if resuming
        progress_file = f".{original_filename}_upload_progress.json"
        progress_path = os.path.join(project.dir, progress_file)
        
        if resume and os.path.exists(progress_path):
            print("Found existing progress file, attempting to resume...")
            # The add_datasource_polars method will handle the resumption logic
        
        # Add the datasource to the project with resume flag
        dodgy_columns = project.add_datasource_polars(
            name=original_filename,
            dataframe=filepath,
            add_to_view=view,
            supplied_columns_only=supplied_only,
            replace_data=replace,
            separator=","
        )
        
        # Get and return the metadata
        metadata = project.get_datasource_metadata(original_filename)
        
        result = {"success": True, "metadata": metadata}
        if dodgy_columns:
            result["dodgy_columns"] = dodgy_columns
            result["warning"] = f"Some columns could not be processed: {dodgy_columns}"
        
        return result
        
    except Exception as e:
        print(f"Error in datasource_processing: {str(e)}")
        raise ValidationError(str(e), status_code=400)

def anndata_processing(project, filepath, original_filename):
    """
    Process an AnnData file (.h5ad) for the project.
    
    Args:
        project (Any): The project instance to which the AnnData will be added.
        filepath (str): Path to the .h5ad file.
        original_filename (str): Original name of the uploaded file.

    Returns:
        dict: A dictionary containing 'success' and message on success.
        
    Raises:
        ValidationError: If any exception occurs during the processing, encapsulating
                         the error message with a 400 status code.
    """
    print(f"Processing AnnData file: {original_filename}")
    print(f"Filepath: {filepath}")

    try:
        # Read the AnnData file
        anndata = sc.read(filepath)
        
        # Convert to MDV format
        convert_scanpy_to_mdv(project.dir, anndata)
        
        result = {"success": True, "message": "AnnData processed successfully"}
        return result
        
    except Exception as e:
        print(f"Error in anndata_processing: {str(e)}")
        raise ValidationError(str(e), status_code=400) from e

def mdv_project_processing(app, projects_base_dir, filepath, original_filename, project_name=None):
    """
    Process a zip file containing an MDV project.
    
    Args:
        app (Flask): Flask application instance to serve the project.
        projects_base_dir (str): Base directory where projects are stored.
        filepath (str): Path to the zip file.
        original_filename (str): Original name of the uploaded file.
        project_name (str, optional): Name for the imported project.

    Returns:
        dict: A dictionary containing 'success', 'project_id', and 'project_name' on success.
        
    Raises:
        ValidationError: If any exception occurs during the processing, encapsulating
                         the error message with a 400 status code.
    """
    import tempfile
    import zipfile
    import shutil
    from mdvtools.mdvproject import MDVProject
    from mdvtools.dbutils.mdv_server_app import is_valid_mdv_project
    from mdvtools.dbutils.dbservice import ProjectService
    
    print(f"Processing MDV project zip file: {original_filename}")
    print(f"Filepath: {filepath}")

    project_path = None
    try:
        # Get next available project ID
        from mdvtools.dbutils.dbservice import ProjectService
        next_id = ProjectService.get_next_project_id()
        if next_id is None:
            raise ValidationError("Failed to determine next project ID from database")

        project_path = os.path.join(projects_base_dir, str(next_id))
        os.makedirs(project_path, exist_ok=True)
        
        # Using a temp directory for extracting files
        with tempfile.TemporaryDirectory() as temp_dir:
            with zipfile.ZipFile(filepath, 'r') as zip_file:
                # Reject entries with absolute paths or ".."
                for file in zip_file.infolist():
                    if file.filename.startswith('/') or '..' in file.filename:
                        raise ValidationError(f"Invalid ZIP file: unsafe path detected - {file.filename}")
                
                temp_extract_path = os.path.join(temp_dir, "extracted")
                os.makedirs(temp_extract_path, exist_ok=True)
                
                # Extract zip file in temp directory
                zip_file.extractall(temp_extract_path)

                extracted_list = os.listdir(temp_extract_path)
                extracted_project_path = None

                # Check if the directory or sub-directory is valid mdv project
                if is_valid_mdv_project(temp_extract_path):
                    extracted_project_path = temp_extract_path
                else:
                    for e in extracted_list:
                        sub_path = os.path.join(temp_extract_path, e)
                        if os.path.isdir(sub_path) and is_valid_mdv_project(sub_path):
                            extracted_project_path = sub_path
                            break

                # Copy the files to newly created project path, if project is valid
                if extracted_project_path is not None:
                    shutil.copytree(extracted_project_path, project_path, dirs_exist_ok=True)
                else:
                    # Clean up the created directory
                    if os.path.exists(project_path):
                        shutil.rmtree(project_path)
                    raise ValidationError("The uploaded file is not a valid MDV project")

        # Create a new MDV project out of the new path and files copied
        mdv_project = MDVProject(project_path, backend_db=True)
        mdv_project.set_editable(True)
        
        # Serve the project through the Flask app
        print(f"Serving new project {next_id} through Flask app")
        mdv_project.serve(app=app, open_browser=False, backend_db=True)
        
        # Initialize the project and register it using project name if valid
        final_project_name = project_name if project_name else os.path.splitext(original_filename)[0]
        new_project = ProjectService.add_new_project(path=project_path, name=final_project_name)

        if not new_project:
            # Clean up on failure
            if os.path.exists(project_path):
                shutil.rmtree(project_path)
            raise ValidationError("Failed to register project in database")
        
        print(f"Successfully created and served project {new_project.id} at /project/{new_project.id}/")
        
        result = {
            "success": True, 
            "message": "MDV project imported successfully",
            "project_id": new_project.id,
            "project_name": new_project.name
        }
        return result
        
    except ValidationError:
        # Re-raise ValidationError as-is
        raise
    except Exception as e:
        print(f"Error in mdv_project_processing: {str(e)}")
        # Clean up project directory if it was created
        if project_path and os.path.exists(project_path):
            try:
                shutil.rmtree(project_path)
            except Exception:
                pass
        raise ValidationError(f"Failed to process MDV project: {str(e)}", status_code=400)
    
class ValidationError(Exception):
    """
    Custom exception raised when datasource validation fails.

    Attributes:
        message (str): Description of the error.
        status_code (int): HTTP status code to associate with the error.
    """
    def __init__(self, message, status_code=400):
        super().__init__(message)
        self.status_code = status_code

def validate_datasource(project, request_data):
    """
    Validates the incoming request data for a datasource upload.

    Ensures that the project is editable, the required fields are present,
    and checks whether a datasource with the same name already exists unless replacement is allowed.

    Args:
        project (Any): Project instance containing state and datasources.
        request_data (dict): Data from the client to validate.

    Returns:
        dict: Validated and normalized request parameters (name, view, replace, supplied_only).

    Raises:
        ValidationError: If the project is read-only or required fields are missing/invalid.
    """
    if "permission" not in project.state or project.state["permission"] != "edit":
        raise ValidationError("Project is read-only")

    name = request_data.get("name")
    if not name:
        raise ValidationError("Request must contain 'name'")

    view = request_data.get("view") if "view" in request_data else "default"

    replace = request_data.get("replace", False)  # Get replace flag (defaults to False)

    if not replace and name in [ds["name"] for ds in project.datasources]:
        raise ValidationError(
            f"Datasource '{name}' already exists, and 'replace' was not set in request"
        )

    supplied_only = request_data.get("supplied_only", False) if "supplied_only" in request_data else False
    
    return {"name": name, "view": view, "replace": replace, "supplied_only": supplied_only}

def validate_anndata(project, request_data):
    """
    Validates the incoming request data for an AnnData upload.

    Ensures that the project is editable.

    Args:
        project (Any): Project instance containing state.
        request_data (dict): Data from the client to validate.

    Returns:
        dict: Validated and normalized request parameters.

    Raises:
        ValidationError: If the project is read-only.
    """
    if "permission" not in project.state or project.state["permission"] != "edit":
        raise ValidationError("Project is read-only")
    
    return {}

def validate_mdv_project(app, projects_base_dir, request_data):
    """
    Validates the incoming request data for an MDV project zip upload.

    Ensures that the projects_base_dir is accessible and validates project name if provided.

    Args:
        app (Flask): Flask application instance.
        projects_base_dir (str): Base directory where projects are stored.
        request_data (dict): Data from the client to validate.

    Returns:
        dict: Validated and normalized request parameters.

    Raises:
        ValidationError: If validation fails.
    """
    if not os.path.exists(projects_base_dir):
        raise ValidationError("Projects base directory not accessible")
    
    # Get project name if provided
    project_name = request_data.get("project_name")
    if project_name:
        # Basic validation for project name
        if len(project_name.strip()) == 0:
            raise ValidationError("Project name cannot be empty")
        # Remove any potentially dangerous characters
        import re
        if not re.match(r'^[a-zA-Z0-9_\-\s]+$', project_name):
            raise ValidationError("Project name contains invalid characters")
    
    return {"project_name": project_name}
