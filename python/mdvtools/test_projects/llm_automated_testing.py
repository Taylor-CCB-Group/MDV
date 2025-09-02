import os
import pandas as pd
import scanpy as sc
import re

import time
import logging
from contextlib import contextmanager
import os
#from openpyxl import load_workbook


from mdvtools.mdvproject import MDVProject

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

from mdvtools.llm.local_files_utils import crawl_local_repo, extract_python_code_from_py, extract_python_code_from_ipynb
from mdvtools.llm.templates import get_createproject_prompt_RAG, prompt_data
from mdvtools.llm.code_manipulation import parse_view_name, prepare_code
from mdvtools.llm.code_execution import execute_code
from mdvtools.llm.chatlog import LangchainLoggingHandler, log_chat_item

import matplotlib

### VARIABLE SETTINGS:
PROJECT_PATH = "mdv/automation/"
DATASET_PATH = "mdv/automation/ilc_viz_ready_revised.h5ad"
QUESTION_LIST_PATH = "python/mdvtools/test_projects/automation_queries_results.xlsx"


# create logger
logger = logging.getLogger('timing_results')
logger.setLevel(logging.DEBUG)

matplotlib.use('Agg') # this should prevent it making any windows etc

print('# setting keys and variables')
# .env file should have OPENAI_API_KEY
load_dotenv()

print('# Crawl the local repository to get a list of relevant file paths')

code_files_urls = crawl_local_repo()

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

embeddings = OpenAIEmbeddings(
    model="text-embedding-3-large"  # Specify the model to use for generating embeddings
    )

db = FAISS.from_documents(texts, embeddings)

retriever = db.as_retriever(
    search_type="similarity",      # Specify the search type as "similarity"
    search_kwargs={"k": 5},        # Set search parameters, in this case, return the top 5 results
)

#***************** Project #*****************

project_path = os.path.expanduser(PROJECT_PATH)

# Load data
data_path = DATASET_PATH
adata = sc.read_h5ad(data_path)
cells_df = pd.DataFrame(adata.obs)
datasource_name = "datasource_name"

# Create project
projectMK = MDVProject(project_path, delete_existing=False)

# Add datasource
projectMK.add_datasource(datasource_name, cells_df)

#***************** LLM #*****************

logger= logging.Logger(__name__)
langchain_logging_handler = LangchainLoggingHandler(logger)
log = logger.info
df1 = None
if len(projectMK.datasources) == 0:
    raise ValueError("The project does not have any datasources")
elif len(projectMK.datasources) > 1:
    ds_name1 = projectMK.datasources[1]['name']
    df1 = projectMK.get_datasource_as_dataframe(ds_name1)
ds_name = projectMK.datasources[0]['name']
try:
    df = projectMK.get_datasource_as_dataframe(ds_name)
    code_llm = ChatOpenAI(temperature=0.1, model="gpt-4o")
    dataframe_llm = ChatOpenAI(temperature=0.1, model="gpt-4o")
    if len(projectMK.datasources) == 1:
        agent = lp.create_pandas_dataframe_agent( # handle_parsing_errors no longer supported
            dataframe_llm, df, verbose=True, allow_dangerous_code=True,
            # handle_parsing_errors="Error in pandas agent"
        )
    elif len(projectMK.datasources) == 2:
        assert(df1 is not None), "The second datasource dataframe is None"
        agent = lp.create_pandas_dataframe_agent(
            dataframe_llm, [df, df1], verbose=True, allow_dangerous_code=True,
            # handle_parsing_errors="Error in pandas agent"
        )
except Exception as e:
    # raise ValueError(f"An error occurred while trying to create the agent: {e[:100]}")
    # todo keep better track of the state of the agent, what went wrong etc
    # pjt - is it ok to continue if the agent creation fails? e.g. hit 'possibly unbound code_llm' below
    raise e
    ok = False

question_file = pd.read_excel(QUESTION_LIST_PATH, sheet_name="Sheet1")
question_list = question_file["queries"].tolist()

pandas_input_list = []
pandas_output_list = []
final_code_list = []
stdout_list = []
stderr_list = []

for question in question_list:
    question = str(question)
    full_prompt = prompt_data + "\nQuestion: " + question
    logger.info(f"Question asked by user: {question}")
    try:
        # Argument of type "str" cannot be assigned to parameter "input" of type "Dict[str, Any]"
        response = agent.invoke(full_prompt) # type: ignore for now
        pandas_input_list.append(response['input'])
        pandas_output_list.append(response['output'])
        #assert('output' in response) # we might allow

        # List all files in the directory
        files_in_dir = os.listdir(projectMK.dir)

        # Initialize variables
        csv_file = None
        h5ad_file = None

        # Identify the CSV or H5AD file
        for file in files_in_dir:
            if file.endswith(".csv"):
                csv_file = file
            elif file.endswith(".h5ad"):
                h5ad_file = file

        # Determine the path_to_data
        if csv_file:
            path_to_data = os.path.join(projectMK.dir, csv_file)
        elif h5ad_file:
            path_to_data = os.path.join(projectMK.dir, h5ad_file)
        else:
            raise FileNotFoundError("No CSV or H5AD file found in the directory.")

        datasource_name = ds_name

        prompt_RAG = get_createproject_prompt_RAG(projectMK, path_to_data, datasource_name, response['output'], question)
        prompt_RAG_template = PromptTemplate(
            template=prompt_RAG,
            input_variables=["context", "question"]
        )

        qa_chain = RetrievalQA.from_llm(
            llm=code_llm,
            prompt=prompt_RAG_template,
            retriever=retriever,
            return_source_documents=True
        )
        context = retriever
        output = qa_chain.invoke({"context": context, "query": question})
        result = output["result"]

        final_code = prepare_code(result, df, projectMK, log, modify_existing_project=True, view_name=question)


        ok, stdout, stderr = execute_code(final_code, open_code=False, log=log)
        final_code_list.append(final_code)
        stdout_list.append(stdout)
        stderr_list.append(stderr)

        if not ok:
            print(f"# Error: code execution failed\n> {stderr}")
        else:
            print(final_code)
            # we need to update this so we don't have the method on the project object
            view_name = parse_view_name(final_code)
            if view_name is None:
                raise Exception("Parsing view name failed")
            log_chat_item(projectMK, question, output, prompt_RAG, final_code, conversation_id="test", context="", view_name=view_name)

            # we want to know the view_name to navigate to as well... for now we do that in the calling code
            #return f"I ran some code for you:\n\n```python\n{final_code}```"
    except Exception as e:
        print(f"Error: {str(e)[:100]}")

pandas_input_list = [
    re.search(r'Question:.*', x).group() if re.search(r'Question.*', x) else "None" # type: ignore - actually, this will raise an error rather than return "None"
    for x in pandas_input_list]

print(question_list)
print(pandas_input_list)
print(pandas_output_list)
print(final_code_list)
print(stdout_list)
print(stderr_list)

# Convert new data into a DataFrame
new_data = pd.DataFrame({"pandas input response": pandas_input_list, "pandas output response": pandas_output_list, "final code": final_code_list, "output stdout": stdout_list, "output stderr": stderr_list})

with pd.ExcelWriter(QUESTION_LIST_PATH, mode='a', if_sheet_exists="overlay") as writer:
    new_data.to_excel(writer, sheet_name='Sheet2')

projectMK.set_editable(True)
projectMK.serve()
