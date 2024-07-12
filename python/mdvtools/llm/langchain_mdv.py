# Code Generation using Retrieval Augmented Generation + LangChain

# importing packages

import os
import pandas as pd
import regex as re
from typing import Optional

from langchain_openai import ChatOpenAI
from langchain_openai import OpenAIEmbeddings
from langchain.schema.document import Document
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import FAISS
from langchain.text_splitter import Language
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate
import langchain_experimental.agents.agent_toolkits.pandas.base as lp
from dotenv import load_dotenv
import subprocess
import tempfile

from .github_utils import crawl_github_repo, extract_python_code_from_py, extract_python_code_from_ipynb

print('# setting keys and variables')
# .env file should have OPENAI_API_KEY & GITHUB_TOKEN
load_dotenv()
# OPENAI_API_KEY environment variable will be used internally by OpenAI modules

mypath = os.path.dirname(__file__)


def extract_code_from_response(response: str):
    """Extracts Python code from a markdown string response."""
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
import json

def load_data(path):
    #Load data from the specified CSV file.
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    #Convert plot data to JSON format.
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))
    
"""

print('# Crawl the GitHub repository to get a list of relevant file URLs')
code_files_urls = crawl_github_repo()

# Initialize an empty list to store the extracted code documents
code_strings = []

# populate code_strings with the code from .py & .ipynb files in code_files_urls
for i in range(0, len(code_files_urls)):
    if code_files_urls[i].endswith(".py"):
        content = extract_python_code_from_py(code_files_urls[i])
        doc = Document(page_content=content, metadata={"url": code_files_urls[i], "file_index": i})
        code_strings.append(doc)
    elif code_files_urls[i].endswith(".ipynb"):
        content_ipynb = extract_python_code_from_ipynb(code_files_urls[i])
        doc_ipynb = Document(page_content=content_ipynb, metadata={"url": code_files_urls[i], "file_index": i})
        code_strings.append(doc_ipynb)

print('# Initialize a text splitter for chunking the code strings')
text_splitter = RecursiveCharacterTextSplitter.from_language(
    language=Language.PYTHON,  # Specify the language as Python
    chunk_size=20000,           # Set the chunk size to 1500 characters
    chunk_overlap=2000          # Set the chunk overlap to 150 characters
)

print('# Split the code documents into chunks using the text splitter')
texts = text_splitter.split_documents(code_strings)

# Set the number of queries per minute (QPM) for embedding requests
# EMBEDDING_QPM = 100 # unused

# Set the number of batches for processing embeddings
# EMBEDDING_NUM_BATCH = 5 # unused

print('# Initialize an instance of the OpenAIEmbeddings class')
embeddings = OpenAIEmbeddings(
    model="text-embedding-3-large"  # Specify the model to use for generating embeddings
    )

print('# Create an index from the embedded code chunks')
print('# Use FAISS (Facebook AI Similarity Search) to create a searchable index')
db = FAISS.from_documents(texts, embeddings)

print('# Initialize the retriever from the FAISS index')
retriever = db.as_retriever(
    search_type="similarity",      # Specify the search type as "similarity"
    search_kwargs={"k": 5},        # Set search parameters, in this case, return the top 5 results
)


def project_wizard(user_question: Optional[str], project_name: str = 'project'):

    print('# Initialize an instance of the ChatOpenAI class with specified parameters')
    # Set the temperature to 0.1 for more deterministic responses
    # Specify the model to use as "gpt-4o"

    code_llm = ChatOpenAI(temperature=0.1, model_name="gpt-4o")

    dataframe_llm = ChatOpenAI(temperature=0.1, model_name="gpt-4")

    #user_question = "Create a heatmap plot of the localisation status vs the UTR length"
    if user_question is None:
        user_question = input("What would you like to ask the LLM?")

    path_to_data = os.path.join(mypath, "sample_data/data_cells.csv")
    df = pd.read_csv(path_to_data)

    df_short = df#.dropna().iloc[:2,1:]

    agent = lp.create_pandas_dataframe_agent(
        dataframe_llm, df_short, verbose=True, handle_parse_errors=True, allow_dangerous_code=True
    )


    full_prompt = prompt + "\nQuestion: " + user_question

    print('# the agent might raise an error. Sometimes repeating the same prompt helps...')
    response = agent.invoke(full_prompt) # agent.run is deprecated
    assert('output' in response)
    final_answer = response['output']

    # todo use a different project foler - nb, this  isn't an f-string at the moment, and has some other template stuff I don't know about
    # what is the difference between ProjectTemplate input_variables and the prompt string itself?
    prompt_RAG = """ You can only generate python code based on the provided context. If a response cannot be formed strictly using the context, politely say you need more information to generate the plot.

    Context: {context}]

    The collection of Python scripts provided in the context, is designed to generate various types of data visualizations 
    using the mdvtools library. Each script focuses on a specific type of plot and follows a common structure that includes loading 
    data from a CSV file, creating a plot using specific parameters, and serving the visualization through an MDV project. 

    All scripts in the context share a common workflow:

    Setup: Define the project path, data path, and view name, the project path should always be: project_path = os.path.expanduser('~/mdv/"""+project_name+"""')
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

    # Define the context for the question (this should be retrieved by the retriever, but showing as an example)
    context = retriever

    print('# Invoke the QA chain with the query and context')
    output = qa_chain.invoke({"context": context, "query": user_question})
    result = output["result"]

    # Extracting the file urls retrieved from the context 

    # unused
    # context_information = output['source_documents']
    # context_information_metadata = [context_information[i].metadata for i in range(len(context_information))]
    # context_information_metadata_url = [context_information_metadata[i]['url'] for i in range(len(context_information_metadata))]
    # context_information_metadata_name = [s[82:] for s in context_information_metadata_url]

    code = extract_code_from_response(result)

    original_script = code
    print('# Apply the reorder transformation')
    modified_script = reorder_parameters(original_script, path_to_data)
    lines = modified_script.splitlines()

    # Find the starting line index
    start_index = next((i for i, line in enumerate(lines) if line.strip().startswith('def')), None)

    if start_index is not None:
        # Capture all lines starting from the first 'def'
        captured_lines = "\n".join(lines[start_index:])
        #print("Captured part:\n", captured_lines)
    else:
        print("Pattern not found")

    # with open("temp_code_3.py", "w") as f:
    #     f.write(packages_functions)
    #     f.write(captured_lines)
        #f.write("\n".join(lines))

    # Log the prompt and the output of the LLM to the google sheets
    # log_to_google_sheet(sheet, str(context_information_metadata_name), output['query'], prompt_RAG, code)

    # print('# Run the saved Python file. This will start a server on localhost:5050, open the browser and display the plot with the server continuing to run in the background.')
    # %run temp_code_3.py
    print('# Executing the code...')
    final_code = f"""{packages_functions}\n{captured_lines}
else:
    main()"""
    final_code = final_code.replace("project.serve()", "# project.serve()")
    print(final_code)
    # exec(final_code) - I don't understand why this doesn't work
    # "name 'MDVProject' is not defined"
    with tempfile.NamedTemporaryFile(suffix='.py') as temp:
        print(f'Writing to temp file: {temp.name}')
        temp.write(final_code.encode())
        temp.flush()
        # subprocess.run(["code", temp.name])
        subprocess.run(["python", temp.name])


if __name__ == "__main__":
    project_wizard(None)