import os
import argparse
import logging
import pandas as pd
import scanpy as sc
import sys
from dotenv import load_dotenv

from mdvtools.mdvproject import MDVProject
from langchain_openai import ChatOpenAI, OpenAIEmbeddings
from langchain.schema.document import Document
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import FAISS
from langchain.text_splitter import Language
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate
import langchain_experimental.agents.agent_toolkits.pandas.base as lp

from mdvtools.llm.local_files_utils import crawl_local_repo, extract_python_code_from_py, extract_python_code_from_ipynb
from mdvtools.llm.templates import get_createproject_prompt_RAG, prompt_data
from mdvtools.llm.code_manipulation import prepare_code
from mdvtools.llm.code_execution import execute_code

import matplotlib
matplotlib.use('Agg')  # Prevent GUI issues


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
    load_dotenv()
    
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
    
    project_path = os.path.expanduser(project_path)
    logger.info(f"Loading dataset from {dataset_path}")
    adata = sc.read_h5ad(dataset_path)
    cells_df = pd.DataFrame(adata.obs)
    
    logger.info("Creating MDV project...")
    project = MDVProject(project_path, delete_existing=False)
    project.add_datasource("datasource_name", cells_df)
    
    logger.info("Setting up LLM interaction...")
    code_llm = ChatOpenAI(temperature=0.1, model="gpt-4o")
    dataframe_llm = ChatOpenAI(temperature=0.1, model="gpt-4o")
    
    ds_name = project.datasources[0]['name']
    df = project.get_datasource_as_dataframe(ds_name)
    agent = lp.create_pandas_dataframe_agent(dataframe_llm, df, verbose=True, allow_dangerous_code=True)
    
    question_file = pd.read_csv(question_list_path, skipinitialspace=True, index_col=False)
    question_list = question_file['requests'].tolist()
    
    results = []
    for question in question_list:
        try:
            response = agent.invoke(prompt_data + "\nQuestion: " + question)
            prompt_RAG = get_createproject_prompt_RAG(project, dataset_path, ds_name, response['output'])
            prompt_template = PromptTemplate(template=prompt_RAG, input_variables=["context", "question"])
            
            qa_chain = RetrievalQA.from_llm(llm=code_llm, prompt=prompt_template, retriever=retriever, return_source_documents=True)
            output = qa_chain.invoke({"context": retriever, "query": question})
            result_code = prepare_code(output["result"], df, project, logger.info, modify_existing_project=True, view_name=question)
            
            ok, stdout, stderr = execute_code(result_code, open_code=False, log=logger.info)
            results.append({
                "question": question,
                "pandas_input": response.get('input', ''),
                "pandas_output": response.get('output', ''),
                "context": response.get('context', ''),
                "final_code": result_code,
                "stdout": stdout,
                "stderr": stderr,
                "execution_status": "Success" if ok else "Failed"
            })
            if not ok:
                logger.error(f"Execution failed: {stderr}")
        except Exception as e:
            logger.error(f"Error processing question: {str(e)}")
    
    result_df = pd.DataFrame(results)
    result_df.to_csv(output_csv, index=False)
    logger.info(f"Results saved to {output_csv}")
    
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run MDV LLM analysis from command line.")
    parser.add_argument("project_path", help="Path to the MDV project directory")
    parser.add_argument("dataset_path", help="Path to the H5AD dataset file")
    parser.add_argument("question_list_path", help="Path to the Excel file containing questions")
    parser.add_argument("output_csv", help="Path to save the output CSV file")
    
    args = parser.parse_args()
    main(args.project_path, args.dataset_path, args.question_list_path, args.output_csv)
