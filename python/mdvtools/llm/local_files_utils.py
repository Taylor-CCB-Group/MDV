import os
import nbformat

# Get the directory of the current script
mypath = os.path.dirname(__file__)

# Define the relative path
DIRECTORY_PATH = os.path.join(mypath, "../test_projects/TAURUS_examples/")

def crawl_local_repo(
    directory_path: str = DIRECTORY_PATH
):
    """
    Crawls a local directory to retrieve file paths based on specified criteria.

    Args:
        directory_path (str): The path to the local project directory.

    Returns:
        list: List of file paths that match the criteria.
    """

    # List of files to ignore
    ignore_list = ["__init__.py", "pbmc3k_tutorial.ipynb", "pbmc3k_tutorial.py"]

    # Initialize an empty list to store file URLs
    files = []


    # Walk through the directory tree
    for root, dirs, file_names in os.walk(os.path.abspath(directory_path)):
        # Skip hidden directories (those starting with '.')
        dirs[:] = [d for d in dirs if not d.startswith('.')]
        
        for file_name in file_names:
            # Check if the file meets the criteria for inclusion
            if file_name not in ignore_list and (file_name.endswith('.py') or file_name.endswith('.ipynb')):
                file_path = os.path.join(root, file_name)
                files.append(file_path)

    # Return the list of collected file paths
    return files

# Extracts the Python code from a .py file from the local filesystem
def extract_python_code_from_py(local_file_path):
    # Read the Python file from the local file system
    with open(local_file_path, 'r', encoding='utf-8') as f:
        python_code = f.read()

    # Return the extracted Python code
    return python_code


# Extracts the Python code from a .ipynb (Jupyter Notebook) file from the local filesystem
def extract_python_code_from_ipynb(local_file_path, cell_type="code"):
    # Read the notebook content from the local file
    with open(local_file_path, 'r', encoding='utf-8') as f:
        notebook_content = f.read()

    # Parse the notebook content using nbformat
    notebook = nbformat.reads(notebook_content, as_version=nbformat.NO_CONVERT)

    # Initialize a variable to store the extracted Python code
    python_code = None

    # Iterate over the cells in the notebook
    for cell in notebook.cells:
        # Check if the cell type matches the specified type
        if cell.cell_type == cell_type:
            # Append the cell's source code to the python_code variable
            if not python_code:
                python_code = cell.source
            else:
                python_code += "\n" + cell.source

    # Return the extracted Python code
    return python_code