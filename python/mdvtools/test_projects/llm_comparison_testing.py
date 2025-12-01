"""
Comparison framework for testing both PythonAstREPLTool and structured output approaches.
This allows for systematic comparison of the two methods.
"""

import os
import pandas as pd
import scanpy as sc
import time
import logging
from contextlib import contextmanager
import traceback
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
from datetime import datetime
import json

from mdvtools.mdvproject import MDVProject
from mdvtools.llm.langchain_mdv import ProjectChat
from mdvtools.llm.structured_project_chat import StructuredProjectChat
from mdvtools.llm.chat_protocol import ChatRequest

from langchain_openai import ChatOpenAI
from langchain_openai import OpenAIEmbeddings
from langchain.schema.document import Document
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import FAISS
from langchain.text_splitter import Language
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate

from dotenv import load_dotenv
from mdvtools.llm.local_files_utils import crawl_local_repo, extract_python_code_from_py, extract_python_code_from_ipynb
from mdvtools.llm.templates import get_createproject_prompt_RAG, prompt_data
from mdvtools.llm.code_manipulation import parse_view_name, prepare_code
from mdvtools.llm.code_execution import execute_code
from mdvtools.llm.chatlog import LangchainLoggingHandler, log_chat_item

import matplotlib

# Configuration
PROJECT_PATH = "mdv/automation/"
DATASET_PATH = "mdv/automation/ilc_viz_ready_revised.h5ad"
QUESTION_LIST_PATH = "python/mdvtools/test_projects/automation_queries_results.xlsx"
OUTPUT_PATH = "python/mdvtools/test_projects/comparison_results.json"

# Setup logging
logger = logging.getLogger('comparison_testing')
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler('comparison_testing.log')
file_handler.setLevel(logging.INFO)
logger.addHandler(file_handler)

matplotlib.use('Agg')

@dataclass
class TestResult:
    """Result of a single test"""
    question: str
    method: str  # "pandas_agent" or "structured"
    success: bool
    execution_time: float
    code_generated: str
    view_name: str
    error_message: Optional[str] = None
    structured_result: Optional[Dict[str, Any]] = None
    pandas_agent_result: Optional[Dict[str, Any]] = None


@dataclass
class ComparisonResult:
    """Result of comparing two methods"""
    question: str
    pandas_agent_result: TestResult
    structured_result: TestResult
    success_match: bool
    execution_time_diff: float
    code_similarity: float  # Placeholder for future implementation


@contextmanager
def time_block(name):
    start_time = time.time()
    yield
    end_time = time.time()
    duration = end_time - start_time
    logger.info(f"Block '{name}' took {duration:.4f} seconds")
    print(f"Block '{name}' took {duration:.4f} seconds")


class ComparisonFramework:
    """Framework for comparing PythonAstREPLTool vs structured output approaches"""
    
    def __init__(self, project_path: str, dataset_path: str):
        self.project_path = project_path
        self.dataset_path = dataset_path
        self.project = None
        self.pandas_chat = None
        self.structured_chat = None
        self.results: List[ComparisonResult] = []
        
        # Initialize RAG components (shared between both methods)
        self._setup_rag()
        
    def _setup_rag(self):
        """Setup RAG components for code generation"""
        print('# Setting up RAG components')
        load_dotenv()
        
        with time_block("RAG setup"):
            # Crawl repository
            code_files_urls = crawl_local_repo()
            code_strings = []
            
            for i, file_path in enumerate(code_files_urls):
                if file_path.endswith(".py"):
                    content = extract_python_code_from_py(file_path)
                    doc = Document(page_content=content, metadata={"url": file_path, "file_index": i})
                    code_strings.append(doc)
                elif file_path.endswith(".ipynb"):
                    content_ipynb = extract_python_code_from_ipynb(file_path)
                    doc_ipynb = Document(page_content=content_ipynb, metadata={"url": file_path, "file_index": i})
                    code_strings.append(doc_ipynb)
            
            # Create embeddings and retriever
            text_splitter = RecursiveCharacterTextSplitter.from_language(
                language=Language.PYTHON,
                chunk_size=20000,
                chunk_overlap=2000
            )
            texts = text_splitter.split_documents(code_strings)
            
            embeddings = OpenAIEmbeddings(model="text-embedding-3-large")
            self.db = FAISS.from_documents(texts, embeddings)
            self.retriever = self.db.as_retriever(
                search_type="similarity",
                search_kwargs={"k": 5},
            )
    
    def setup_project(self):
        """Setup the MDV project and both chat instances"""
        print('# Setting up project and chat instances')
        
        with time_block("Project setup"):
            # Load data
            adata = sc.read_h5ad(self.dataset_path)
            cells_df = pd.DataFrame(adata.obs)
            datasource_name = "datasource_name"
            
            # Create project
            self.project = MDVProject(self.project_path, delete_existing=False)
            self.project.add_datasource(datasource_name, cells_df)
            
            # Create both chat instances
            self.pandas_chat = ProjectChat(self.project)
            self.structured_chat = StructuredProjectChat(self.project)
    
    def run_single_test(self, question: str, method: str) -> TestResult:
        """Run a single test with the specified method"""
        print(f'Running {method} test for: {question}')
        
        start_time = time.time()
        
        try:
            # Create a mock chat request
            chat_request = ChatRequest(
                id="test_id",
                room="test_room",
                conversation_id="test_conversation",
                message=question,
                handle_error=lambda error, extra_metadata=None: print(f"Error: {error}")
            )
            assert self.pandas_chat is not None
            assert self.structured_chat is not None
            # Run the appropriate method
            if method == "pandas_agent":
                result = self.pandas_chat.ask_question(chat_request)
                assert result is not None
                structured_result = None
                pandas_agent_result = result
            else:  # structured
                result = self.structured_chat.ask_question(chat_request)
                assert result is not None
                structured_result = result.get("structured_result")
                pandas_agent_result = None
            
            execution_time = time.time() - start_time
            assert pandas_agent_result is not None
            assert structured_result is not None
            return TestResult(
                question=question,
                method=method,
                success=not result.get("error", True),
                execution_time=execution_time,
                code_generated=result.get("code", "") or "",
                view_name=result.get("view_name", "") or "",
                error_message=result.get("message") if result.get("error") else None,
                structured_result=structured_result,
                pandas_agent_result=pandas_agent_result # type: ignore
            )
            
        except Exception as e:
            execution_time = time.time() - start_time
            error_message = f"{str(e)[:500]}\n\n{traceback.format_exc()}"
            
            return TestResult(
                question=question,
                method=method,
                success=False,
                execution_time=execution_time,
                code_generated="",
                view_name="",
                error_message=error_message
            )
    
    def compare_methods(self, question: str) -> ComparisonResult:
        """Compare both methods for a single question"""
        print(f'Comparing methods for: {question}')
        
        # Run pandas agent test
        pandas_result = self.run_single_test(question, "pandas_agent")
        
        # Run structured test
        structured_result = self.run_single_test(question, "structured")
        
        # Calculate comparison metrics
        success_match = pandas_result.success == structured_result.success
        execution_time_diff = structured_result.execution_time - pandas_result.execution_time
        
        # Simple code similarity (could be enhanced with more sophisticated comparison)
        code_similarity = 0.0
        if pandas_result.code_generated and structured_result.code_generated:
            # Basic similarity based on common imports and structure
            pandas_code = pandas_result.code_generated.lower()
            structured_code = structured_result.code_generated.lower()
            
            # Count common elements
            common_elements = 0
            total_elements = 0
            
            # Check for common imports
            imports = ["import pandas", "import scanpy", "from mdvtools", "def main"]
            for imp in imports:
                if imp in pandas_code and imp in structured_code:
                    common_elements += 1
                total_elements += 1
            
            if total_elements > 0:
                code_similarity = common_elements / total_elements
        
        return ComparisonResult(
            question=question,
            pandas_agent_result=pandas_result,
            structured_result=structured_result,
            success_match=success_match,
            execution_time_diff=execution_time_diff,
            code_similarity=code_similarity
        )
    
    def run_comparison_suite(self, questions: List[str]) -> List[ComparisonResult]:
        """Run comparison tests for a list of questions"""
        print(f'Running comparison suite for {len(questions)} questions')
        
        results = []
        
        for i, question in enumerate(questions):
            print(f'Processing question {i+1}/{len(questions)}: {question}')
            
            try:
                comparison_result = self.compare_methods(question)
                results.append(comparison_result)
                
                # Log results
                logger.info(f"Question {i+1}: {question}")
                logger.info(f"  Pandas Agent: Success={comparison_result.pandas_agent_result.success}, "
                          f"Time={comparison_result.pandas_agent_result.execution_time:.2f}s")
                logger.info(f"  Structured: Success={comparison_result.structured_result.success}, "
                          f"Time={comparison_result.structured_result.execution_time:.2f}s")
                logger.info(f"  Success Match: {comparison_result.success_match}")
                logger.info(f"  Time Difference: {comparison_result.execution_time_diff:.2f}s")
                logger.info(f"  Code Similarity: {comparison_result.code_similarity:.2f}")
                
            except Exception as e:
                logger.error(f"Error processing question {i+1}: {e}")
                print(f"Error processing question {i+1}: {e}")
                continue
        
        self.results = results
        return results
    
    def generate_summary_report(self) -> Dict[str, Any]:
        """Generate a summary report of the comparison results"""
        if not self.results:
            return {"error": "No results to summarize"}
        
        total_questions = len(self.results)
        pandas_successes = sum(1 for r in self.results if r.pandas_agent_result.success)
        structured_successes = sum(1 for r in self.results if r.structured_result.success)
        
        pandas_avg_time = sum(r.pandas_agent_result.execution_time for r in self.results) / total_questions
        structured_avg_time = sum(r.structured_result.execution_time for r in self.results) / total_questions
        
        success_matches = sum(1 for r in self.results if r.success_match)
        avg_code_similarity = sum(r.code_similarity for r in self.results) / total_questions
        
        return {
            "summary": {
                "total_questions": total_questions,
                "pandas_agent_success_rate": pandas_successes / total_questions,
                "structured_success_rate": structured_successes / total_questions,
                "pandas_agent_avg_time": pandas_avg_time,
                "structured_avg_time": structured_avg_time,
                "success_match_rate": success_matches / total_questions,
                "avg_code_similarity": avg_code_similarity
            },
            "detailed_results": [
                {
                    "question": r.question,
                    "pandas_agent": {
                        "success": r.pandas_agent_result.success,
                        "execution_time": r.pandas_agent_result.execution_time,
                        "view_name": r.pandas_agent_result.view_name,
                        "error": r.pandas_agent_result.error_message
                    },
                    "structured": {
                        "success": r.structured_result.success,
                        "execution_time": r.structured_result.execution_time,
                        "view_name": r.structured_result.view_name,
                        "error": r.structured_result.error_message
                    },
                    "comparison": {
                        "success_match": r.success_match,
                        "execution_time_diff": r.execution_time_diff,
                        "code_similarity": r.code_similarity
                    }
                }
                for r in self.results
            ]
        }
    
    def save_results(self, output_path: str):
        """Save results to a JSON file"""
        report = self.generate_summary_report()
        report["timestamp"] = datetime.now().isoformat()
        
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"Results saved to {output_path}")


def main():
    """Main function to run the comparison framework"""
    print("Starting LLM Method Comparison Framework")
    
    # Load questions
    try:
        question_file = pd.read_excel(QUESTION_LIST_PATH, sheet_name="Sheet1")
        question_list = question_file["queries"].tolist()
        print(f"Loaded {len(question_list)} questions from {QUESTION_LIST_PATH}")
    except Exception as e:
        print(f"Error loading questions: {e}")
        # Fallback to a few test questions
        question_list = [
            "Can you create a scatter plot?",
            "Show me a histogram",
            "Create a bar chart"
        ]
        print(f"Using fallback questions: {question_list}")
    
    # Initialize framework
    framework = ComparisonFramework(PROJECT_PATH, DATASET_PATH)
    
    try:
        # Setup project
        framework.setup_project()
        
        # Run comparison suite
        results = framework.run_comparison_suite(question_list)
        
        # Generate and save report
        framework.save_results(OUTPUT_PATH)
        
        # Print summary
        report = framework.generate_summary_report()
        summary = report["summary"]
        
        print("\n" + "="*50)
        print("COMPARISON SUMMARY")
        print("="*50)
        print(f"Total Questions: {summary['total_questions']}")
        print(f"Pandas Agent Success Rate: {summary['pandas_agent_success_rate']:.2%}")
        print(f"Structured Success Rate: {summary['structured_success_rate']:.2%}")
        print(f"Pandas Agent Avg Time: {summary['pandas_agent_avg_time']:.2f}s")
        print(f"Structured Avg Time: {summary['structured_avg_time']:.2f}s")
        print(f"Success Match Rate: {summary['success_match_rate']:.2%}")
        print(f"Avg Code Similarity: {summary['avg_code_similarity']:.2%}")
        print("="*50)
        
    except Exception as e:
        print(f"Error in comparison framework: {e}")
        logger.error(f"Framework error: {e}")
        traceback.print_exc()


if __name__ == "__main__":
    main() 