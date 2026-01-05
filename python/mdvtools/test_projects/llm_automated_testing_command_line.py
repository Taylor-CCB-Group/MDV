import os
import argparse
import logging
import pandas as pd
import scanpy as sc
import sys
from dotenv import load_dotenv
import re
import tempfile
import shutil

from mdvtools.mdvproject import MDVProject
from mdvtools.conversions import convert_scanpy_to_mdv
from langchain_openai import ChatOpenAI, OpenAIEmbeddings
from langchain.schema.document import Document
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import FAISS
from langchain.text_splitter import Language
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate

# packages for custom langchain agent
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain_experimental.tools.python.tool import PythonAstREPLTool
from langchain.agents import create_openai_functions_agent, AgentExecutor

# packages for history
from langchain.memory import ConversationBufferMemory
from langchain.agents import initialize_agent

from mdvtools.llm.local_files_utils import extract_python_code_from_py, extract_python_code_from_ipynb
from mdvtools.llm.templates import get_createproject_prompt_RAG, prompt_data
from mdvtools.llm.code_manipulation import prepare_code
from mdvtools.llm.code_execution import execute_code

from langchain.chains import LLMChain
from pydantic.v1 import BaseModel, Field


import matplotlib
matplotlib.use('Agg')  # Prevent GUI issues

# ******************************************************
# Get the directory of the current script
mypath = os.path.dirname(__file__)

# Define the relative path
DIRECTORY_PATH = os.path.join(mypath, "../test_projects/RAG_examples/ANNDATA_examples/")

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
    ignore_list = ["__init__.py"]

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
# ******************************************************

def setup_logging(log_file="mdv_llm.log"):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),  # Save logs to a file
            logging.StreamHandler(sys.stdout)  # Also print logs to console
        ]
    )
    return logging.getLogger(__name__)


def main(project_path, dataset_path, question_list_path, output_csv):
    logger = setup_logging()
    load_dotenv(os.path.join(os.path.dirname(__file__), "../llm/.env"))

    
    logger.info("Crawling the local repository...")
    code_files_urls = crawl_local_repo()
    
    code_strings = []
    for file_path in code_files_urls:
        if file_path.endswith(".py"):
            content = extract_python_code_from_py(file_path)
        elif file_path.endswith(".ipynb"):
            content = extract_python_code_from_ipynb(file_path)
        else:
            continue
        doc = Document(page_content=content, metadata={"url": file_path})
        code_strings.append(doc)
    
    logger.info("Splitting code documents into chunks...")
    text_splitter = RecursiveCharacterTextSplitter.from_language(
        language=Language.PYTHON, chunk_size=20000, chunk_overlap=2000
    )
    texts = text_splitter.split_documents(code_strings)
    
    embeddings = OpenAIEmbeddings(model="text-embedding-3-large")
    db = FAISS.from_documents(texts, embeddings)
    retriever = db.as_retriever(search_type="similarity", search_kwargs={"k": 5})
    
    temp_dir = None
    if project_path is None or dataset_path is None:
        temp_dir = tempfile.mkdtemp()
        project_path = temp_dir
        logger.info("Using temporary directory for project path.")

        # Load sample data
        logger.info("Loading sample dataset pbmc3k.")
        adata = sc.datasets.pbmc3k()

        logger.info("Creating MDV project...")
        project = convert_scanpy_to_mdv(project_path, adata)
    else:
        project_path = os.path.expanduser(project_path)
        logger.info(f"Loading dataset from {dataset_path}")
        adata = sc.read_h5ad(dataset_path)

        logger.info("Creating MDV project...")
        project = convert_scanpy_to_mdv(project_path, adata)

    logger.info("Setting up LLM interaction...")
    code_llm = ChatOpenAI(temperature=0.1, model="gpt-4o")
    dataframe_llm = ChatOpenAI(temperature=0.1, model="gpt-4o")
    
    if len(project.datasources) >= 2:
        df_list = [project.get_datasource_as_dataframe(ds['name']) for ds in project.datasources[:2]]
        logger.info("Using two data sources for the pandas agent.")
    else:
        df_list = [project.get_datasource_as_dataframe(project.datasources[0]['name'])]
        logger.info("Using one data source for the pandas agent.")

    datasource_names = [ds['name'] for ds in project.datasources[:2]]  # Get names of up to 2 datasources

    #CUSTOM AGENT:
    def create_custom_pandas_agent(llm, dfs: dict, prompt_data, verbose=False):
        """
        Creates a LangChain agent that can interact with Pandas DataFrames using a Python REPL tool.
        
        :param llm: The LLM to use (e.g., OpenAI GPT-4).
        :param dfs: A dictionary of named Pandas DataFrames.
        :param verbose: If True, prints debug information.
        :return: An agent that can answer questions about the DataFrames.
        """
        
        # Step 1: Initialize Memory with Chat History
        memory = ConversationBufferMemory(memory_key="chat_history", return_messages=True)
        
        # Step 2: Create the Python REPL Tool
        python_tool = PythonAstREPLTool(locals={}, globals={})#, sanitize_input=True)

        assert(python_tool.globals is not None)
        python_tool.globals.update(dfs)

        python_tool.globals["list_globals"] = lambda: list(python_tool.globals.keys()) # type: ignore

        # Step 3: Define Contextualization Chain

        contextualize_q_system_prompt = """Given a chat history and the latest user question \
        which might reference context in the chat history, formulate a standalone question \
        which can be understood without the chat history. Do NOT answer the question, \
        just reformulate it if needed and otherwise return it as is. \
        """

        contextualize_prompt = ChatPromptTemplate.from_messages([
            ("system", contextualize_q_system_prompt),
            ("human", "Chat History:\n{chat_history}\n\nUser Question:\n{input}"),])
        
        contextualize_chain = LLMChain(llm=llm, prompt=contextualize_prompt, memory=memory)

        # Step 4: Define the Agent Prompt

        prompt_data_template = f"""You have access to the following Pandas DataFrames: 
        {', '.join(dfs.keys())}. These are preloaded, so do not redefine them.
        Before answering any user question, you must first run `df1.columns` and `df2.columns` to inspect available fields. 
        Use these to correct the column names mentioned by the user.
        You must always invoke the PythonAstREPLTool to check the DataFrames columns and explore the values of the DataFrames.
        Use `df.info()` or `df.index()`. 
        Before running any code, check available variables using `list_globals()`.""" + prompt_data

        prompt = ChatPromptTemplate.from_messages([
            ("system", prompt_data_template),
            MessagesPlaceholder(variable_name="chat_history"),
            ("human", prompt_data_template + "{input}"),
            ("ai", "{agent_scratchpad}"),
        ])

        # Step 5: Create the Agent
        agent = create_openai_functions_agent(llm, [python_tool], prompt)

        # Step 6: Wrap in an Agent Executor (Finalized Agent)
        agent_executor = AgentExecutor(agent=agent, tools=[python_tool], memory=memory, verbose=verbose, return_intermediate_steps=True)

        # Step 7: Wrapper Function to Use Contextualization and Preserve Memory
        def agent_with_contextualization(question):
            #logger.info(f"\n===== New Question =====\n{question}")
            standalone_question = contextualize_chain.run(input=question)
            logger.info(f"Reformulated question: {standalone_question}")
            # Point 1: Log reformulation
            # Point 2: Log what you're sending to the agent
            response = agent_executor.invoke({"input": standalone_question})
            logger.info(f"Agent raw response: {response}")
            # Point 3: Log agent output
            memory.save_context({"input": question}, {"output": response.get("output", str(response))})
            return response
        
        # Attach the memory object to the agent so it can be cleared externally
        agent_with_contextualization.memory = memory # type: ignore not sure why pyright complains about this

        return agent_with_contextualization

    question_file = pd.read_csv(question_list_path, skipinitialspace=True, index_col=False)
    question_list = question_file['requests'].tolist()
    results = []

    dataframes_for_agent = {"df1": pd.DataFrame(df_list[0]), "df2": pd.DataFrame(df_list[1])}

    agent = create_custom_pandas_agent(dataframe_llm, dataframes_for_agent, prompt_data, verbose=True)

    RESET_THRESHOLD = 1  # Set the desired number of requests after which memory should be reset.
    request_counter = 0

    for question in question_list:
        try:
            request_counter += 1  # Increase request counter at each iteration.
            full_prompt = prompt_data + "\nQuestion: " + question
            #logger.info(f"Agent prompt input:\n{full_prompt}")

            ## custom agent code

            response = agent(full_prompt)
            logger.info(f"Pandas agent output: {response.get('output', '')}")

            match = re.search(r'charts\s+(.*)', response['output'])
            charts_part = match.group(1) if match else response['output']

            
            prompt_RAG = get_createproject_prompt_RAG(project, dataset_path, datasource_names[0], response['output'], response['input'])
            #logger.info(f"RAG prompt being sent:\n{prompt_RAG}")

            prompt_template = PromptTemplate(template=prompt_RAG, input_variables=["context", "question"])
            
            qa_chain = RetrievalQA.from_llm(llm=code_llm, prompt=prompt_template, retriever=retriever, return_source_documents=True)
            output = qa_chain.invoke({"query": charts_part})#({"query": response['input'] + response['output']}) #"context": retriever,
            logger.info(f"Raw output:\n{output}")

            result_code = prepare_code(output["result"], df_list[0], project, logger.info, modify_existing_project=True, view_name=question)
            logger.info(f"Generated code:\n{result_code}")

            context_information = output['source_documents']
            context_information_metadata = [context_information[i].metadata for i in range(len(context_information))]
            context_information_metadata_url = [context_information_metadata[i]['url'] for i in range(len(context_information_metadata))]
            
            ok, stdout, stderr = execute_code(result_code, open_code=False, log=logger.info)

            results.append({
                "question": question,
                "pandas_input": response.get('prompt_data', ''), #'input', ''),
                "pandas_output": response.get('output', ''),
                "RAG_output": output,
                "context": context_information_metadata_url,
                "final_code": result_code,
                "stdout": stdout,
                "stderr": stderr,
                "execution_status": "Success" if ok else "Failed"
            })
            if not ok:
                logger.error(f"Execution failed: {stderr}")
        except Exception as e:
            logger.error(f"Error processing question: {str(e)}")

        # Reset the agent's conversation memory after RESET_THRESHOLD requests.
        if request_counter % RESET_THRESHOLD == 0:
            agent.memory.clear()  # Clear the conversation history. # type: ignore
            logger.info(f"Agent memory has been reset after {request_counter} requests.")
    
    result_df = pd.DataFrame(results)
    result_df.to_csv(output_csv, index=False)
    logger.info(f"Results saved to {output_csv}")
    
    project.set_editable(True)
    project.serve()

    # Clean up the temporary directory if it was used
    if temp_dir is not None and project_path == temp_dir:
        shutil.rmtree(temp_dir)
        logger.info("Temporary directory cleaned up.")

if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Run MDV LLM analysis from command line.")
    # parser.add_argument("project_path", help="Path to the MDV project directory")
    # parser.add_argument("dataset_path", help="Path to the H5AD dataset file")
    # parser.add_argument("question_list_path", help="Path to the Excel file containing questions")
    # parser.add_argument("output_csv", help="Path to save the output CSV file")
    
    # args =parser.parse_args()
    #main(args.project_path, args.dataset_path, args.question_list_path, args.output_csv)
    # todo - have a default version of the questions file in the repo.
    main(None, None, f"chat_testing/logs/pbmc3k_questions.csv", "chat_testing/logs/pbmc3k_544_latest.csv")







