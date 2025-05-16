def datasource_processing(project, filepath, original_filename, view, replace, supplied_only):
    """
    Processes a CSV datasource by adding it to the project and retrieving metadata.

    Args:
        project (Any): The project instance to which the datasource will be added.
        filepath (str): Path to the CSV file.
        original_filename (str): Original name of the uploaded file.
        view (str): The view to associate the datasource with.
        replace (bool): Whether to replace existing data.
        supplied_only (bool): Whether to use only supplied columns (currently unused in call).

    Returns:
        dict: Dictionary containing 'success' and metadata or error message.
    """
    print(f"Processing datasource: {original_filename} in view: {view}")
    print(f"Filepath: {filepath}")
    print(f"Replace: {replace}, Supplied only: {supplied_only}")

    try:
        # Add the datasource to the project
        project.add_datasource(
            name=original_filename,
            dataframe=filepath,
            add_to_view=view,
            supplied_columns_only=supplied_only,
            replace_data=replace,
            separator=","
        )
        # Get and return the metadata
        metadata = project.get_datasource_metadata(original_filename)
        return {"success": True, "metadata": metadata}
        
    except Exception as e:
        print(f"Error in datasource_processing: {str(e)}")
        import traceback
        traceback.print_exc()
        return {"success": False, "error": str(e)}, 400
    
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
