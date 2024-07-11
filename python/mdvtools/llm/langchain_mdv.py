# %% [markdown]
# # Code Generation using Retrieval Augmented Generation + LangChain

# %%
# importing packages

import nbformat
import os
import pandas as pd
import requests
import time
import regex as re
from typing import Optional

from langchain_openai import ChatOpenAI
from langchain_openai import OpenAIEmbeddings
from langchain.schema.document import Document
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.vectorstores import FAISS
from langchain.text_splitter import Language
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate
import langchain_experimental.agents.agent_toolkits.pandas.base as lp
from dotenv import load_dotenv

# %%
print('# setting keys and variables')
# .env file should have OPENAI_API_KEY & GITHUB_TOKEN
# also currently need a key.json file with google-sheet config
load_dotenv()
# OPENAI_API_KEY environment variable will be used internally by OpenAI modules
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN")

GITHUB_REPO = "Taylor-CCB-Group/MDV" # @param {type:"string"}
BRANCH_NAME = "mk-API"
PROJECT_PATH_1 = "python/mdvtools/charts"
PROJECT_PATH_2 = "python/mdvtools/test_projects"

FILE_URL_PATH = "./"
FILE_URL_NAME = "code_files_URL.txt"

mypath = os.path.dirname(__file__)

def crawl_github_repo(url: str, is_sub_dir: bool, branch_name: str, project_path: str, access_token=f"{GITHUB_TOKEN}"):
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
    ignore_list = ['__init__.py', 'pbmc3k_tutorial.ipynb', 'pbmc3k_tutorial.py']

    # Determine the appropriate API URL based on whether it's a sub-directory
    if not is_sub_dir:
        api_url = f"https://api.github.com/repos/{url}/contents/{project_path}?ref={branch_name}"
    else:
        api_url = url

    # Set up headers for the GitHub API request, including authorization
    headers = {
        "Accept": "application/vnd.github.v3+json",
        "Authorization": f"Bearer {access_token}"
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
        if item['type'] == 'file' and item['name'] not in ignore_list and (item['name'].endswith('.py') or item['name'].endswith('.ipynb')):
            files.append(item['html_url'])
        # Check if the item is a directory (excluding hidden ones)
        elif item['type'] == 'dir' and not item['name'].startswith("."):
            # Recursively crawl the sub-directory
            sub_files = crawl_github_repo(item['url'], True, branch_name, project_path)
            # Pause briefly to avoid rate limiting
            time.sleep(.1)
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

def extract_code_from_response(response: str):
    """Extracts Python code from a string response."""
    # Use a regex pattern to match content between triple backticks
    code_pattern = r"```python(.*?)```"
    match = re.search(code_pattern, response, re.DOTALL)

    if match:
        # Extract the matched code and strip any leading/trailing whitespaces
        return match.group(1).strip()
    return None

def reorder_parameters(script: str, dataframe_path: str):
    # Load the dataframe to infer the column types
    df = pd.read_csv(dataframe_path)
    categorical_columns = df.select_dtypes(
        include=["object", "category"]
    ).columns.tolist()
    numerical_columns = df.select_dtypes(include=["number"]).columns.tolist()

    def is_categorical(column):
        return column in categorical_columns

    def is_numerical(column):
        return column in numerical_columns

    # Define a regex pattern to find function definitions that create BoxPlots
    patterns = [
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)BoxPlot\((.*?)\)", re.DOTALL
        ),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)DotPlot\((.*?)\)", re.DOTALL
        ),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)AbundanceBoxPlot\((.*?)\)",
            re.DOTALL,
        ),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)HistogramPlot\((.*?)\)", re.DOTALL
        ),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)RingChart\((.*?)\)", re.DOTALL
        ),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)RowChart\((.*?)\)", re.DOTALL
        ),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)StackedRowChart\((.*?)\)", re.DOTALL
        ),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)HeatmapPlot\((.*?)\)", re.DOTALL
        ),
    ]

    for pattern in patterns:
        if pattern.search(script):
            # Define a regex pattern to find params and param patterns
            pattern_param = re.compile(r'params\s*=\s*\[.*?\]|param\s*=\s*".*?"')

            def reorder_params(match_param):
                matched_text = match_param.group(0)  # Get the entire matched text

                # Extract parameter names
                if "params" in matched_text:
                    param_list = re.findall(r"\'(.*?)\'", matched_text)
                else:
                    param_list = [re.findall(r"\'(.*?)\'", matched_text)[0]]

                # Check for the presence of categorical and numerical variables
                has_categorical = any(is_categorical(param) for param in param_list)
                has_numerical = any(is_numerical(param) for param in param_list)

                # Add a categorical variable if none is present
                if not has_categorical and categorical_columns:
                    param_list.insert(0, categorical_columns[0])
                    has_categorical = True

                if len(param_list) < 2:
                    return matched_text  # No need to reorder if there are fewer than 2 parameters

                first_param = param_list[0]
                second_param = param_list[1]

                # Check the types of the parameters using the dataframe
                # if first_param in df.columns and second_param in df.columns:
                if has_categorical and has_numerical:
                    if not (
                        is_categorical(first_param) and is_numerical(second_param)
                    ):
                        param_list[0], param_list[1] = param_list[1], param_list[0]

                # Reconstruct the parameters with reordered values
                if "params" in matched_text:
                    reordered_params = (
                        f"params = ['{param_list[0]}', '{param_list[1:]}']"
                    )
                else:
                    reordered_params = f'param = "{param_list[0]}"'

                return reordered_params.replace("'[", " ").replace("]'", "")

            # Substitute the matches with reordered parameters
            modified_script = re.sub(pattern_param, reorder_params, script)

            return modified_script

    return script

prompt = """
Based on the question asked and the CSV, please perform the following steps:

1. Identify the data asked for in the question.
2. Based on step 1, find the relevant column names in the CSV based on the information identified earlier in the question asked regarding data.
3. Provide the relevant column names from step 2 in a list.
4. For the purposes of this query, if the plot required is one of these:
    a. heatmap
    b. dot plot
    c. box plot
the suitable parameters should always be a categorigal variable (object or categorical) followed by numerical variables.
"""


packages_functions = """import os
import pandas as pd
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.heatmap_plot import HeatmapPlot
from mdvtools.charts.histogram_plot import HistogramPlot
from mdvtools.charts.dot_plot import DotPlot
from mdvtools.charts.box_plot import BoxPlot
from mdvtools.charts.scatter_plot_3D import ScatterPlot3D
from mdvtools.charts.row_chart import RowChart
from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.abundance_box_plot import AbundanceBoxPlot
from mdvtools.charts.stacked_row_plot import StackedRowChart
from mdvtools.charts.ring_chart import RingChart
from mdvtools.charts.violin_plot import ViolinPlot

ViolinPlot

import json \n
\n

def load_data(path):
    #Load data from the specified CSV file.
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    #Convert plot data to JSON format.
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))
    
"""



def project_wizard(user_question: Optional[str]):

    # %%
    print('# Initialize an instance of the ChatOpenAI class with specified parameters')
    # Set the temperature to 0.1 for more deterministic responses
    # Specify the model to use as "gpt-4o"

    code_llm = ChatOpenAI(temperature=0.1, model_name="gpt-4o")

    dataframe_llm = ChatOpenAI(temperature=0.1, model_name="gpt-4")



    # %%
    print('# Crawl the GitHub repository to get a list of relevant file URLs')
    code_files_urls = crawl_github_repo(GITHUB_REPO, False, BRANCH_NAME, PROJECT_PATH_2, GITHUB_TOKEN)

    # Write the list of file URLs to a specified text file
    with open(FILE_URL_PATH + FILE_URL_NAME, 'w') as f:
        # Iterate through the list of file URLs and write each one to the file
        for item in code_files_urls:
            f.write(item + '\n')

    # %%


    # Read the list of file URLs from the specified text file
    with open(FILE_URL_PATH + FILE_URL_NAME) as f:
        code_files_urls = f.read().splitlines()

    # %%
    # Initialize an empty list to store the extracted code documents
    code_strings = []

    # Iterate over the list of file URLs
    for i in range(0, len(code_files_urls)):
        # Check if the file URL ends with ".py"
        if code_files_urls[i].endswith(".py"):
            # Extract the Python code from the .py file
            content = extract_python_code_from_py(code_files_urls[i])
            # Create a Document object with the extracted content and metadata
            doc = Document(page_content=content, metadata={"url": code_files_urls[i], "file_index": i})
            # Append the Document object to the code_strings list
            code_strings.append(doc)
            # Check if the file URL ends with ".py"
        elif code_files_urls[i].endswith(".ipynb"):
            # Extract the Python code from the .py file
            content_ipynb = extract_python_code_from_ipynb(code_files_urls[i])
            # Create a Document object with the extracted content and metadata
            doc_ipynb = Document(page_content=content_ipynb, metadata={"url": code_files_urls[i], "file_index": i})
            # Append the Document object to the code_strings list
            code_strings.append(doc_ipynb)

    # %%
    print('# Initialize a text splitter for chunking the code strings')
    text_splitter = RecursiveCharacterTextSplitter.from_language(
        language=Language.PYTHON,  # Specify the language as Python
        chunk_size=20000,           # Set the chunk size to 1500 characters
        chunk_overlap=2000          # Set the chunk overlap to 150 characters
    )

    print('# Split the code documents into chunks using the text splitter')
    texts = text_splitter.split_documents(code_strings)

    # %%
    # Set the number of queries per minute (QPM) for embedding requests
    EMBEDDING_QPM = 100

    # Set the number of batches for processing embeddings
    EMBEDDING_NUM_BATCH = 5

    print('# Initialize an instance of the OpenAIEmbeddings class')
    embeddings = OpenAIEmbeddings(
        model="text-embedding-3-large"  # Specify the model to use for generating embeddings
        )

    # %%
    print('# Create an index from the embedded code chunks')
    print('# Use FAISS (Facebook AI Similarity Search) to create a searchable index')
    db = FAISS.from_documents(texts, embeddings)

    # %%
    print('# Initialize the retriever from the FAISS index')
    retriever = db.as_retriever(
        search_type="similarity",      # Specify the search type as "similarity"
        search_kwargs={"k": 5},        # Set search parameters, in this case, return the top 5 results
    )

    # %%
    #user_question = "Create a heatmap plot of the localisation status vs the UTR length"
    if user_question is None:
        user_question = input("What would you like to ask the LLM?")

    # %%
    path_to_data = os.path.join(mypath, "sample_data/data_cells.csv")
    df = pd.read_csv(path_to_data)

    # %%
    df_short = df#.dropna().iloc[:2,1:]

    # %%
    agent = lp.create_pandas_dataframe_agent(
        dataframe_llm, df_short, verbose=True, handle_parse_errors=True, allow_dangerous_code=True
    )


    full_prompt = prompt + "\nQuestion: " + user_question

    print('# the agent might raise an error. Sometimes repeating the same prompt helps...')
    final_answer = agent.run(full_prompt)

    # %%
    prompt_RAG = """ You can only generate python code based on the provided context. If a response cannot be formed strictly using the context, politely say you need more information to generate the plot.

    Context: {context}]

    The collection of Python scripts provided in the context, is designed to generate various types of data visualizations 
    using the mdvtools library. Each script focuses on a specific type of plot and follows a common structure that includes loading 
    data from a CSV file, creating a plot using specific parameters, and serving the visualization through an MDV project. 

    All scripts in the context share a common workflow:

    Setup: Define the project path, data path, and view name, the project path should always be: project_path = os.path.expanduser('~/mdv/project')
    Plot function definition: Define the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) using a function in the same way as the context.
    Project Creation: Initialize an MDVProject instance using the method: MDVProject(project_path, delete_existing=True).
    Data Loading: Load data from the specified CSV file into a pandas DataFrame using the load_data(path) function.
    Data adding: Add the data source to the project using the method: project.add_datasource(data_path, data).
    Plot Creation: Create the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) and define the plot paramaters in the same way as in the context.
    Data Conversion: Convert the plot data to JSON format for integration with the MDV project using the convert_plot_to_json(plot) function.
    Serving: Configure the project view, set it to editable, and serve the project using the .set_view(view_name, plot_view), .set_editable(True) and .serve() methods.

    You are a top-class Python developer. Generate a Python script following the workflow detailed above and use exactly the same lines of code as the scripts in the context. 
    The data should be loaded from a csv provided, the path to the csv is given by: """ + path_to_data + """ 
    This list """ + final_answer + """ specifies the names of the data fields that need to be plotted, for example in the params field. Get the structure of params definition from the context. 
    The question: {question} specifies the plot required. 
    """

    #The plot you should create is the same as the plot created in the context. Specify the parameters according to the respective files in the context for each plot type. DO NOT add any parameters that have not been defined previously.

    # %%
    print('# Create a PromptTemplate object using the defined RAG prompt')
    prompt_RAG_template = PromptTemplate(
        template=prompt_RAG,          # Specify the template string
        input_variables=["context", "question"]  # Define the input variables for the template
    )

    # Initialize a RetrievalQA chain using the specified language model, prompt template, and retriever
    qa_chain = RetrievalQA.from_llm(
        llm=code_llm,                 # Specify the language model to use
        prompt=prompt_RAG_template,   # Use the defined prompt template
        retriever=retriever,          # Use the initialized retriever for context retrieval
        return_source_documents=True  # Configure the chain to return source documents
    )

    # %%
    # Define the context for the question (this should be retrieved by the retriever, but showing as an example)
    context = retriever

    print('# Invoke the QA chain with the query and context')
    output = qa_chain.invoke({"context": context, "query": user_question})
    result = output["result"]

    # %%
    # Extracting the file urls retrieved from the context 

    context_information = output['source_documents']
    context_information_metadata = [context_information[i].metadata for i in range(len(context_information))]
    context_information_metadata_url = [context_information_metadata[i]['url'] for i in range(len(context_information_metadata))]
    context_information_metadata_name = [s[82:] for s in context_information_metadata_url]

    # %%
    
    code = extract_code_from_response(result)

    # %%
    original_script = code

    # %%

    # %%
    print('# Apply the reorder transformation')
    modified_script = reorder_parameters(original_script, path_to_data)

    # %%


    lines = modified_script.splitlines()

    # Find the starting line index
    start_index = next((i for i, line in enumerate(lines) if line.strip().startswith('def')), None)

    if start_index is not None:
        # Capture all lines starting from the first 'def'
        captured_lines = "\n".join(lines[start_index:])
        #print("Captured part:\n", captured_lines)
    else:
        print("Pattern not found")

    with open("temp_code_3.py", "w") as f:
        f.write(packages_functions)
        f.write(captured_lines)
        #f.write("\n".join(lines))

    # %%
    # Log the prompt and the output of the LLM to the google sheets
    # log_to_google_sheet(sheet, str(context_information_metadata_name), output['query'], prompt_RAG, code)

    # %%
    print('# Run the saved Python file. This will start a server on localhost:5050, open the browser and display the plot with the server continuing to run in the background.')
    # %run temp_code_3.py


if __name__ == "__main__":
    project_wizard(None)