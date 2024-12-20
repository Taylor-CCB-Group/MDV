import time
import logging
from contextlib import contextmanager

# Code Generation using Retrieval Augmented Generation + LangChain
from typing import Callable
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

from .local_files_utils import crawl_local_repo, extract_python_code_from_py, extract_python_code_from_ipynb
from .templates import get_createproject_prompt_RAG, prompt_data, get_updateproject_prompt_RAG
from .code_manipulation import prepare_code
from .code_execution import execute_code

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

print('# Initialize the retriever from the FAISS index')

with time_block("b5: Retriever creating"):
    retriever = db.as_retriever(
        search_type="similarity",      # Specify the search type as "similarity"
        search_kwargs={"k": 5},        # Set search parameters, in this case, return the top 5 results
    )

# ... for testing basic mechanics without invoking the agent
mock_agent = False

class ProjectChat():
    def __init__(self, project: MDVProject, log: Callable[[str], None] = print):
        self.project = project
        self.log = log
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
                    self.agent = lp.create_pandas_dataframe_agent( # handle_parse_errors no longer supported
                        self.dataframe_llm, self.df, verbose=True, allow_dangerous_code=True
                    )
                elif len(project.datasources) == 2:
                    self.agent = lp.create_pandas_dataframe_agent(
                        self.dataframe_llm, [self.df, self.df1], verbose=True, allow_dangerous_code=True
                    )
            self.ok = True
        except Exception as e:
            # raise ValueError(f"An error occurred while trying to create the agent: {e[:100]}")
            log(f"An error occurred while trying to create the agent: {str(e)[:100]}")
            self.ok = False

    def ask_question(self, question: str): # async?
        if mock_agent:
            ok, strdout, stderr = execute_code('import mdvtools\nprint("mdvtools import ok")')
            return f"mdvtools import ok? {ok}\n\nstdout: {strdout}\n\nstderr: {stderr}"
        with time_block("b9a: Pandas agent invoking"):
            # shortly after this...
            # Error in StdOutCallbackHandler.on_chain_start callback: AttributeError("'NoneType' object has no attribute 'get'")
            self.log(f"Asking the LLM: '{question}'")
            if not self.ok:
                return "This agent is not ready to answer questions"
            full_prompt = prompt_data + "\nQuestion: " + question
            print(full_prompt)
            logger.info(f"Question asked by user: {question}")
        try:
            with time_block("b9b: Pandas agent invoking"):
                # Argument of type "str" cannot be assigned to parameter "input" of type "Dict[str, Any]"
                response = self.agent.invoke(full_prompt) # type: ignore for now
                assert('output' in response) # we might allow
            #!!! csv_path is not wanted - the code tries to use that as data source name which is all wrong
            with time_block("b10: RAG prompt preparation"):
                #! todo - allow for the ability to pass a path to some data that may not already be in the project?
                # at least we should allow for multiple datasources. Do we even need to pass ds_name(s), or should the script look it (them) up for itself?
                path_to_data = self.project.dir + "/anndata.h5ad"
                #!!!!!! for now, assuming there will be an anndata.h5ad file in the project directory and will fail ungacefully if there isn't!!!!
                prompt_RAG = get_createproject_prompt_RAG(self.project.id, path_to_data, response['output']) #self.ds_name, response['output'])
                prompt_RAG_template = PromptTemplate(
                    template=prompt_RAG,
                    input_variables=["context", "question"]
                )

            with time_block("b11: RAG chain"):
                qa_chain = RetrievalQA.from_llm(
                    llm=self.code_llm,
                    prompt=prompt_RAG_template,
                    retriever=retriever,
                    return_source_documents=True
                )
                context = retriever
                output = qa_chain.invoke({"context": context, "query": question})
                result = output["result"]

            with time_block("b12: Prepare code"):
                final_code = prepare_code(result, self.df, self.log, modify_existing_project=True, view_name=question)
            # log_to_google_sheet(sheet, str(context_information_metadata_name), output['query'], prompt_RAG, code)
            # todo - save code at various stages of processing...
            # log_chat(output, prompt_RAG, final_code)
            with time_block("b13: Chat logging by MDV"):
                self.project.log_chat_item(output, prompt_RAG, final_code)
            with time_block("b14: Execute code"):
                ok, stdout, stderr = execute_code(final_code, open_code=False, log=self.log)
            if not ok:
                return f"# Error: code execution failed\n> {stderr}"
            else:
                self.log(final_code)
                return f"I ran some code for you:\n\n```python\n{final_code}```"
        except Exception as e:
            return f"Error: {str(e)[:100]}"
