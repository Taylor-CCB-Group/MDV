import time
import logging
from contextlib import contextmanager

# Code Generation using Retrieval Augmented Generation + LangChain
import os
import pandas as pd
from typing import Optional, Callable
from mdvtools.mdvproject import MDVProject

from langchain_openai import ChatOpenAI
from langchain_openai import OpenAIEmbeddings
from langchain.schema.document import Document
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import FAISS
from langchain.text_splitter import Language
from langchain.chains import RetrievalQA # why can't I see the type for this?
from langchain.prompts import PromptTemplate
import langchain_experimental.agents.agent_toolkits.pandas.base as lp
from dotenv import load_dotenv

from .github_utils import crawl_github_repo, extract_python_code_from_py, extract_python_code_from_ipynb
from .templates import get_createproject_prompt_RAG, prompt_data
from .code_manipulation import prepare_code
from .code_execution import execute_code

#logging.basicConfig(filename="timing_results.log", level=logging.INFO)

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
import matplotlib
matplotlib.use('Agg') # this should prevent it making any windows etc

print('# setting keys and variables')
# .env file should have OPENAI_API_KEY & GITHUB_TOKEN
load_dotenv()
# OPENAI_API_KEY environment variable will be used internally by OpenAI modules

mypath = os.path.dirname(__file__)


print('# Crawl the GitHub repository to get a list of relevant file URLs')
with time_block("b1: GitHub repo crawling"):
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

with time_block("b2: Text splitter initialising"):
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
with time_block("b3: Embeddings creating"):
    embeddings = OpenAIEmbeddings(
        model="text-embedding-3-large"  # Specify the model to use for generating embeddings
        )

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

class ProjectChat():
    def __init__(self, project: MDVProject, log: Callable[[str], None] = print):
        self.project = project
        self.log = log
        if len(project.datasources) == 0:
            raise ValueError("The project does not have any datasources")
        elif len(project.datasources) > 1:
            log("The project has more than one datasource, only the first one will be used")
        self.ds_name = project.datasources[0]['name']
        try:
            self.df = project.get_datasource_as_dataframe(self.ds_name)
            with time_block("b6: Initialising LLM for RAG"):
                self.code_llm = ChatOpenAI(temperature=0.1, model_name="gpt-4o")
            with time_block("b7: Initialising LLM for agent"):
                self.dataframe_llm = ChatOpenAI(temperature=0.1, model_name="gpt-4")
            with time_block("b8: Pandas agent creating"):
                self.agent = lp.create_pandas_dataframe_agent(
                    self.dataframe_llm, self.df, verbose=True, handle_parse_errors=True, allow_dangerous_code=True
                )
            self.ok = True
        except Exception as e:
            # raise ValueError(f"An error occurred while trying to create the agent: {e[:100]}")
            log(f"An error occurred while trying to create the agent: {str(e)[:100]}")
            self.ok = False
    
    def ask_question(self, question: str): # async?
        with time_block("b9a: Pandas agent invoking"):
            self.log(f"Asking the LLM: '{question}'")
            if not self.ok:
                return "This agent is not ready to answer questions"
            full_prompt = prompt_data + "\nQuestion: " + question
            print(full_prompt)
            logger.info(f"Question asked by user: {question}")
        try:
            with time_block("b9b: Pandas agent invoking"):
                response = self.agent.invoke(full_prompt) 
                assert('output' in response)
            #!!! csv_path is not wanted - the code tries to use that as data source name which is all wrong
            with time_block("b10: RAG prompt preparation"):
                prompt_RAG = get_createproject_prompt_RAG(self.project.id, self.ds_name, response['output'])
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
            self.project.log_chat_item(output, prompt_RAG, final_code)
            with time_block("b13: Execute code"):
                execute_code(final_code, open_code=False, log=self.log)
                self.log(final_code)
                return f"I ran some code for you:\n\n```python\n{final_code}```"
        except Exception as e:
            return f"Error: {str(e)[:100]}"


def project_wizard(user_question: Optional[str], project_name: str = 'project', log: Callable[[str], None] = print):

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

    # could we get this to log with the log function?
    agent = lp.create_pandas_dataframe_agent(
        dataframe_llm, df_short, verbose=True, handle_parse_errors=True, allow_dangerous_code=True
    )

    full_prompt = prompt_data + "\nQuestion: " + user_question

    log('# the agent might raise an error. Sometimes repeating the same prompt helps...')
    response = agent.invoke(full_prompt) # agent.run is deprecated
    assert('output' in response)
    final_answer = response['output']

    prompt_RAG = get_createproject_prompt_RAG(project_name, path_to_data, final_answer)

    #The plot you should create is the same as the plot created in the context. Specify the parameters according to the respective files in the context for each plot type. DO NOT add any parameters that have not been defined previously.

    log('# Create a PromptTemplate object using the defined RAG prompt')
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

    log('# Invoke the QA chain with the query and context')
    output = qa_chain.invoke({"context": context, "query": user_question})
    result = output["result"]

    # extract and process the code from the response
    # nb - what was path_to_data is currently only used within reorder_parameters here.
    # the actual string is not used in the generated code, so we can pass df instead
    final_code = prepare_code(result, df, log)
    log(final_code)
    # passing `log` around is a chore, would be nice to know if there's a better way
    execute_code(final_code, log=log)


if __name__ == "__main__":
    project_wizard(None)