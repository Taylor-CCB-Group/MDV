from dotenv import load_dotenv
import requests
import os
import time
import nbformat

load_dotenv()

GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN")

GITHUB_REPO = "Taylor-CCB-Group/MDV"  # @param {type:"string"}
# BRANCH_NAME = "mk-API"
COMMIT_HASH = "af4192b" # latest commit by mk as of this writing
PROJECT_PATH_1 = "python/mdvtools/charts"
PROJECT_PATH_2 = "python/mdvtools/test_projects"


def crawl_github_repo(
    url: str = GITHUB_REPO,
    is_sub_dir: bool = False,
    branch_or_commit_name: str = COMMIT_HASH,
    project_path: str = PROJECT_PATH_2,
    access_token=f"{GITHUB_TOKEN}",
):
    """
    Crawls a GitHub repository to retrieve file URLs based on specified criteria.

    Args:
        url (str): The GitHub repository URL or sub-directory URL.
        is_sub_dir (bool): Flag indicating if the current URL is a sub-directory.
        branch_name (str): The branch name to crawl.
        project_path (str): The path of the project in the repository.
        access_token (str, optional): GitHub access token for authentication. Defaults to GITHUB_TOKEN.

    Returns:
        list: List of file URLs that match the criteria.
    """

    # List of files to ignore
    ignore_list = ["__init__.py", "pbmc3k_tutorial.ipynb", "pbmc3k_tutorial.py"]

    # Determine the appropriate API URL based on whether it's a sub-directory
    if not is_sub_dir:
        api_url = f"https://api.github.com/repos/{url}/contents/{project_path}?ref={branch_or_commit_name}"
    else:
        api_url = url

    # Set up headers for the GitHub API request, including authorization
    headers = {
        "Accept": "application/vnd.github.v3+json",
        "Authorization": f"Bearer {access_token}",
    }

    # Make a GET request to the GitHub API
    response = requests.get(api_url, headers=headers)
    # Raise an exception for any request errors
    response.raise_for_status()

    # Initialize an empty list to store file URLs
    files = []

    # Parse the JSON response content
    contents = response.json()

    # Iterate over the items in the contents
    for item in contents:
        # Check if the item is a file and meets the criteria for inclusion
        if (
            item["type"] == "file"
            and item["name"] not in ignore_list
            and (item["name"].endswith(".py") or item["name"].endswith(".ipynb"))
        ):
            files.append(item["html_url"])
        # Check if the item is a directory (excluding hidden ones)
        elif item["type"] == "dir" and not item["name"].startswith("."):
            # Recursively crawl the sub-directory
            sub_files = crawl_github_repo(item["url"], True, branch_or_commit_name, project_path)
            # Pause briefly to avoid rate limiting
            time.sleep(0.1)
            # Add the sub-directory files to the list
            files.extend(sub_files)

    # Return the list of collected file URLs
    return files


# Extracts the Python code from a .ipynb (Jupyter Notebook) file from GitHub
def extract_python_code_from_ipynb(github_url: str, cell_type="code"):
    # Convert the GitHub URL to the raw content URL
    raw_url = github_url.replace("github.com", "raw.githubusercontent.com").replace(
        "/blob/", "/"
    )

    # Make a GET request to fetch the raw content of the notebook
    response = requests.get(raw_url)
    response.raise_for_status()  # Check for any request errors

    # Get the content of the notebook as text
    notebook_content = response.text

    # Read the notebook content using nbformat
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


# Extracts the Python code from a .py file from GitHub
def extract_python_code_from_py(github_url):
    # Convert the GitHub URL to the raw content URL
    raw_url = github_url.replace("github.com", "raw.githubusercontent.com").replace(
        "/blob/", "/"
    )

    # Make a GET request to fetch the raw content of the Python file
    response = requests.get(raw_url)
    response.raise_for_status()  # Check for any request errors

    # Get the content of the Python file as text
    python_code = response.text
    # print(python_code)

    # Return the extracted Python code
    return python_code
