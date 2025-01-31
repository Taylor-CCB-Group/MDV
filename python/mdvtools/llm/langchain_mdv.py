import time
import logging
from contextlib import contextmanager
import os

# Code Generation using Retrieval Augmented Generation + LangChain
# from typing import Callable
from mdvtools.mdvproject import MDVProject

from langchain_openai import ChatOpenAI
from langchain_openai import OpenAIEmbeddings
from langchain.schema.document import Document
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import FAISS
from langchain.text_splitter import Language
# https://python.langchain.com/api_reference/langchain/chains/langchain.chains.retrieval_qa.base.RetrievalQA.html
# Deprecated since version 0.1.17: This class is deprecated.
# Use the create_retrieval_chain constructor instead.
# See migration guide here: https://python.langchain.com/docs/versions/migrating_chains/retrieval_qa/
from langchain.chains import RetrievalQA # why can't I see the type for this? see above comment...
# from langchain_core.output_parsers import StrOutputParser
# from langchain_core.runnables import RunnablePassthrough

from langchain.prompts import PromptTemplate
import langchain_experimental.agents.agent_toolkits.pandas.base as lp
from dotenv import load_dotenv
from mdvtools.websocket import ChatSocketAPI

from .local_files_utils import crawl_local_repo, extract_python_code_from_py, extract_python_code_from_ipynb
from .templates import get_createproject_prompt_RAG, prompt_data, get_updateproject_prompt_RAG
from .code_manipulation import parse_view_name, prepare_code
from .code_execution import execute_code
from .chatlog import LangchainLoggingHandler

import matplotlib

# create logger
logger = logging.getLogger('timing_results')
logger.setLevel(logging.DEBUG)

# create file handler and set level to INFO
file_handler = logging.FileHandler('timing_results.log')
file_handler.setLevel(logging.INFO)
logger.addHandler(file_handler)
logger.info(time.asctime())

@contextmanager
def time_block(name):
    start_time = time.time()  # Setup: Start timing
    yield  # This is where the code inside the `with` block runs
    end_time = time.time()  # Teardown: After `with` block is done
    duration = end_time - start_time
    # Only log the message if it contains 'Block'
    logger.info(f"Block '{name}' took {duration:.4f} seconds")
    print(f"Block '{name}' took {duration:.4f} seconds")


# sometimes this tries to draw with matplotlib... which causes an error that we can't catch
# we don't control the code that's being run by the agent...
# but we can call `matplotlib.use('Agg')` before calling the agent (and any other code that might try to draw)
matplotlib.use('Agg') # this should prevent it making any windows etc

print('# setting keys and variables')
# .env file should have OPENAI_API_KEY
load_dotenv()
# OPENAI_API_KEY environment variable will be used internally by OpenAI modules

# todo - add some testing so that we notice if we introduce breaking changes like this
# mypath = os.path.dirname(__file__)
# path_to_data = os.path.join(mypath, "sample_data/bcell_viz_ready_revised.h5ad")

print('# Crawl the local repository to get a list of relevant file paths')
with time_block("b1: Local repo crawling"):
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

with time_block("b2: Text splitter initialising"):
    text_splitter = RecursiveCharacterTextSplitter.from_language(
        language=Language.PYTHON,  # Specify the language as Python
        chunk_size=20000,           # Set the chunk size to 1500 characters
        chunk_overlap=2000          # Set the chunk overlap to 150 characters
    )
    print('# Split the code documents into chunks using the text splitter')
    texts = text_splitter.split_documents(code_strings)

#assert(len(texts) > 0)

# Set the number of queries per minute (QPM) for embedding requests
# EMBEDDING_QPM = 100 # unused

# Set the number of batches for processing embeddings
# EMBEDDING_NUM_BATCH = 5 # unused

print('# Initialize an instance of the OpenAIEmbeddings class')
with time_block("b3: Embeddings creating"):
    embeddings = OpenAIEmbeddings(
        model="text-embedding-3-large"  # Specify the model to use for generating embeddings
    )

print(embeddings)
print('# Create an index from the embedded code chunks')
print('# Use FAISS (Facebook AI Similarity Search) to create a searchable index')


with time_block("b4: FAISS database creating"):
    db = FAISS.from_documents(texts, embeddings)

print("# Initialize the retriever from the FAISS index")

with time_block("b5: Retriever creating"):
    retriever = db.as_retriever(
        search_type="similarity",      # Specify the search type as "similarity"
        search_kwargs={"k": 5},        # Set search parameters, in this case, return the top 5 results
    )

# ... for testing basic mechanics without invoking the agent
mock_agent = False


class ProjectChat:
    def __init__(self, project: MDVProject):
        self.project = project
        self.socket_api = ChatSocketAPI(project)
        logger = self.socket_api.logger
        self.langchain_logging_handler = LangchainLoggingHandler(logger)
        self.config = {"callbacks": [self.langchain_logging_handler]}
        log = self.log = logger.info
        if len(project.datasources) == 0:
            raise ValueError("The project does not have any datasources")
        elif len(project.datasources) > 1:
            # tempting to think we should just give the agent all of the datasources here
            log("The project has more than one datasource, only the first one will be used")
            self.ds_name1 = project.datasources[1]['name']
            self.df1 = project.get_datasource_as_dataframe(self.ds_name1)
        self.ds_name = project.datasources[0]['name']
        try:
            self.df = project.get_datasource_as_dataframe(self.ds_name)
            with time_block("b6: Initialising LLM for RAG"):
                self.code_llm = ChatOpenAI(temperature=0.1, model="gpt-4o")
            with time_block("b7: Initialising LLM for agent"):
                self.dataframe_llm = ChatOpenAI(temperature=0.1, model="gpt-4o")
            with time_block("b8: Pandas agent creating"):
                if len(project.datasources) == 1:
                    self.agent = lp.create_pandas_dataframe_agent( # handle_parsing_errors no longer supported
                        self.dataframe_llm, self.df, verbose=True, allow_dangerous_code=True,
                        # handle_parsing_errors="Error in pandas agent"
                    )
                elif len(project.datasources) == 2:
                    self.agent = lp.create_pandas_dataframe_agent(
                        self.dataframe_llm, [self.df, self.df1], verbose=True, allow_dangerous_code=True,
                        # handle_parsing_errors="Error in pandas agent"
                    )
            self.ok = True
        except Exception as e:
            # raise ValueError(f"An error occurred while trying to create the agent: {e[:100]}")
            log(f"An error occurred while trying to create the agent: {str(e)[:100]}")
            # todo keep better track of the state of the agent, what went wrong etc
            self.ok = False

    def ask_question(self, question: str, id: str):  # async?
        """
        Ask a question, generate code to answer it, execute the code...

        How should we stream updates on the progress of the code execution back to the user?
        We have a log, connected to a websocket... but perhaps it would make sense to yield the output of the code execution
        and feed that back to the http response as it comes in?

        Then in the front-end, we can have a view of the verbose logging stuff, but not necessarily
        """
        progress = 0
        if mock_agent:
            ok, strdout, stderr = execute_code(
                'import mdvtools\nprint("mdvtools import ok")'
            )
            return f"mdvtools import ok? {ok}\n\nstdout: {strdout}\n\nstderr: {stderr}"
        with time_block("b9a: Pandas agent invoking"):  # ~0.005% of time
            self.socket_api.update_chat_progress(
                "Pandas agent invoking", id, progress, 0
            )
            # shortly after this...
            # Error in StdOutCallbackHandler.on_chain_start callback: AttributeError("'NoneType' object has no attribute 'get'")
            self.log(f"Asking the LLM: '{question}'")
            if not self.ok:
                return "This agent is not ready to answer questions"
            full_prompt = prompt_data + "\nQuestion: " + question
            print(full_prompt)
            # provide an update to the user that we're working on it, using provided id.
            logger.info(f"Question asked by user: {question}")
        try:
            with time_block("b9b: Pandas agent invoking"):  # ~31.4% of time
                self.socket_api.update_chat_progress(
                    "Pandas agent invoking...", id, progress, 31
                )
                progress += 31
                # Argument of type "str" cannot be assigned to parameter "input" of type "Dict[str, Any]"
                response = self.agent.invoke(
                    full_prompt, config={"callbacks": [self.langchain_logging_handler]}
                )  # type: ignore for now
                assert "output" in response  # we might allow
            with time_block("b10: RAG prompt preparation"):  # ~0.003% of time
                self.socket_api.update_chat_progress(
                    "RAG prompt preparation...", id, progress, 1
                )
                # List all files in the directory
                files_in_dir = os.listdir(self.project.dir)

                # Initialize variables
                csv_file = None
                h5ad_file = None

                # Identify the CSV or H5AD file
                # subject to review at a later date
                for file in files_in_dir:
                    if file.endswith(".csv"):
                        csv_file = file
                    elif file.endswith(".h5ad"):
                        h5ad_file = file

                # Determine the path_to_data
                if csv_file:
                    path_to_data = os.path.join(self.project.dir, csv_file)
                elif h5ad_file:
                    path_to_data = os.path.join(self.project.dir, h5ad_file)
                else:
                    raise FileNotFoundError("No CSV or H5AD file found in the directory.")

                # path_to_data now contains the correct file path
                #path_to_data = "data_cells.csv"
                datasource_name = self.ds_name

                #!!!!!! for now, assuming there will be an anndata.h5ad file in the project directory and will fail ungacefully if there isn't!!!!
                # we pass a reference to the actual project object and let figuring out the path be an internal implementation detail...
                # this should be more robust, and also more flexible in terms of what reasoning this method may be able to do internally in the future
                # I appear to have an issue though - the configuration of the devcontainer doesn't flag whether or not the thing we're passing is the right type
                # and the assert in the function is being triggered even though it should be fine
                prompt_RAG = get_createproject_prompt_RAG(self.project, path_to_data, datasource_name, response['output']) #self.ds_name, response['output'])
                prompt_RAG_template = PromptTemplate(
                    template=prompt_RAG,
                    input_variables=["context", "question"]
                )

            with time_block("b11: RAG chain"):  # ~60% of time
                self.socket_api.update_chat_progress(
                    "Invoking RAG chain...", id, progress, 60
                )
                progress += 60
                qa_chain = RetrievalQA.from_llm(
                    llm=self.code_llm,
                    prompt=prompt_RAG_template,
                    retriever=retriever,
                    return_source_documents=True,
                )
                context = retriever
                output = qa_chain.invoke({"context": context, "query": question})
                result = output["result"]
                print(result)

            with time_block("b12: Prepare code"):  # <0.1% of time
                final_code = prepare_code(
                    result,
                    self.df,
                    self.project,
                    self.log,
                    modify_existing_project=True,
                    view_name=question,
                )
                # view_name = parse_view_name(final_code)
            # log_to_google_sheet(sheet, str(context_information_metadata_name), output['query'], prompt_RAG, code)
            # todo - save code at various stages of processing...
            # log_chat(output, prompt_RAG, final_code)
            with time_block("b13: Chat logging by MDV"):  # <0.1% of time
                self.project.log_chat_item(output, prompt_RAG, final_code)
            with time_block("b14: Execute code"):  # ~9% of time
                self.socket_api.update_chat_progress(
                    "Executing code...", id, progress, 9
                )
                progress += 9
                ok, stdout, stderr = execute_code(
                    final_code, open_code=False, log=self.log
                )
            if not ok:
                return f"# Error: code execution failed\n> {stderr}"
            else:
                self.log(final_code)
                self.socket_api.update_chat_progress(
                    "Finished processing query", id, 100, 0
                )
                # we want to know the view_name to navigate to as well... for now we do that in the calling code
                return f"I ran some code for you:\n\n```python\n{final_code}```"
        except Exception as e:
            return f"Error: {str(e)[:100]}"
