"""
Structured version of ProjectChat that uses structured output instead of PythonAstREPLTool.
This version generates structured data analysis and chart configurations directly.
"""

import time
import logging
from contextlib import contextmanager
import os
import traceback
from typing import List, Dict

from mdvtools.llm.chat_protocol import AskQuestionResult, ChatRequest, ProjectChatProtocol
from mdvtools.mdvproject import MDVProject
from mdvtools.llm.structured_schemas import (
    DataAnalysisResult, ChartGenerationResult, StructuredQueryResult,
    create_chart_config, ChartConfig
)

from langchain_openai import ChatOpenAI
from langchain_openai import OpenAIEmbeddings
from langchain.schema.document import Document
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import FAISS
from langchain.text_splitter import Language
from langchain_core.pydantic_v1 import BaseModel, Field
from langchain.schema import HumanMessage
from langchain.chains import LLMChain
from langchain.memory import ConversationBufferMemory
from langchain.prompts import PromptTemplate
from langchain.chains import RetrievalQA

from dotenv import load_dotenv
from mdvtools.llm.chatlog import ChatSocketAPI, LangchainLoggingHandler, log_chat_item
from mdvtools.llm.local_files_utils import crawl_local_repo, extract_python_code_from_py, extract_python_code_from_ipynb
from mdvtools.llm.templates import get_createproject_prompt_RAG, prompt_data
from mdvtools.llm.code_manipulation import parse_view_name, prepare_code
from mdvtools.llm.code_execution import execute_code
from mdvtools.llm.markdown_utils import create_suggested_questions_prompt

import matplotlib

# create logger
logger = logging.getLogger('timing_results')
logger.setLevel(logging.DEBUG)

# create file handler and set level to INFO
file_handler = logging.FileHandler('timing_results.log')
file_handler.setLevel(logging.INFO)
logger.addHandler(file_handler)
logger.info(time.asctime())

# create a separate logger for chat debugging
chat_debug_logger = logging.getLogger('chat_debug')
chat_debug_logger.setLevel(logging.DEBUG)
chat_debug_handler = logging.FileHandler('chat_debug.log')
chat_debug_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
chat_debug_handler.setFormatter(formatter)
chat_debug_logger.addHandler(chat_debug_handler)


class SuggestedQuestions(BaseModel):
    """A list of suggested questions to ask about the data."""
    questions: List[str] = Field(
        description="A list of 5 suggested questions based on the provided data context. Each question should be a single, complete sentence."
    )


@contextmanager
def time_block(name):
    start_time = time.time()
    yield
    end_time = time.time()
    duration = end_time - start_time
    logger.info(f"Block '{name}' took {duration:.4f} seconds")
    print(f"Block '{name}' took {duration:.4f} seconds")


# Prevent matplotlib from creating windows
matplotlib.use('Agg')

print('# setting keys and variables')
load_dotenv()
if not os.getenv("OPENAI_API_KEY"):
    print("OPENAI_API_KEY is not set in the environment variables. Please set it if LLM integration is required.")
    raise ValueError("OPENAI_API_KEY is not set in the environment variables. Please set it if LLM integration is required.")

print('# Crawl the local repository to get a list of relevant file paths')
with time_block("b1: Local repo crawling"):
    code_files_urls = crawl_local_repo()
    code_strings = []
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
        language=Language.PYTHON,
        chunk_size=20000,
        chunk_overlap=2000
    )
    print('# Split the code documents into chunks using the text splitter')
    texts = text_splitter.split_documents(code_strings)

print('# Initialize an instance of the OpenAIEmbeddings class')
with time_block("b3: Embeddings creating"):
    embeddings = OpenAIEmbeddings(model="text-embedding-3-large")

print('# Create an index from the embedded code chunks')
with time_block("b4: FAISS database creating"):
    db = FAISS.from_documents(texts, embeddings)

print("# Initialize the retriever from the FAISS index")
with time_block("b5: Retriever creating"):
    retriever = db.as_retriever(
        search_type="similarity",
        search_kwargs={"k": 5},
    )


class StructuredProjectChat(ProjectChatProtocol):
    def __init__(self, project: MDVProject):
        self.project = project
        self.welcome = (
            "Hello, I'm an AI assistant that has access to the data in this project "
            "and is designed to help build views for visualising it using structured analysis. What can I help you with?"
        )
        self.log = logger.info
        self.conversation_memories = {}
        self.suggested_questions: list[str] = []
        
        if len(project.datasources) == 0:
            raise ValueError("The project does not have any datasources")
        elif len(project.datasources) > 1:
            self.log("The project has more than one datasource, only the first one will be used")

        self.ds_name = project.datasources[0]['name']
        try:
            try:
                with time_block("b_suggest: Generating suggested questions"):
                    llm = ChatOpenAI(temperature=0.1, model="gpt-4.1")
                    structured_llm = llm.with_structured_output(SuggestedQuestions)
                    prompt_text = create_suggested_questions_prompt(self.project)
                    messages = [HumanMessage(content=prompt_text)]
                    result = structured_llm.invoke(messages)
                    self.suggested_questions = result.questions  # type: ignore
                    self.log(f"Generated suggested questions: {self.suggested_questions}")
            except Exception as e:
                self.log(f"Could not generate suggested questions: {e}")
                self.suggested_questions = [str(e)]

            self.init_error = False
        except Exception as e:
            error_message = f"{str(e)[:500]}\n\n{traceback.format_exc()}"
            self.log(error_message)
            self.error_message = error_message
            self.init_error = True

    def get_suggested_questions(self):
        if self.init_error:
            return []
        if self.suggested_questions:
            return self.suggested_questions
        return ["This should be unreachable"]
    
    def get_or_create_memory(self, conversation_id: str):
        """Get or create conversation memory for a specific conversation"""
        if conversation_id not in self.conversation_memories:
            self.conversation_memories[conversation_id] = ConversationBufferMemory(
                memory_key="chat_history", 
                return_messages=True
            )
        return self.conversation_memories[conversation_id]

    def clear_conversation(self, conversation_id: str):
        """Clear the chat history for a specific conversation"""
        if conversation_id in self.conversation_memories:
            del self.conversation_memories[conversation_id]
        chat_debug_logger.info(f"Cleared conversation history for conversation {conversation_id}")

    def analyze_data_structured(self, llm: ChatOpenAI, question: str, project: MDVProject) -> DataAnalysisResult:
        """Analyze data using structured output based on project metadata instead of DataFrames"""
        
        # Create structured LLM for data analysis
        structured_llm = llm.with_structured_output(DataAnalysisResult)
        
        # Get project metadata in markdown format (same as used for suggested questions)
        from mdvtools.llm.markdown_utils import create_project_markdown
        project_info = create_project_markdown(project)
        
        # Create analysis prompt using project metadata
        analysis_prompt = f"""
        You are analyzing data to help create visualizations. You have access to the following project metadata:
        
        {project_info}
        
        Based on the user question: "{question}"
        
        Please analyze the data and select appropriate columns for visualization. Consider:
        1. What type of visualization is being requested
        2. What columns contain the relevant data
        3. Whether the data is categorical, numerical, or text
        4. How many columns are needed for the requested chart type
        5. Which datasource(s) contain the relevant data
        
        For gene-related queries, look for gene names in the data and use appropriate gene expression columns.
        For categorical data, select columns that represent categories or groups.
        For numerical data, select columns that represent measurements or values.
        
        IMPORTANT: The project may have multiple datasources. A given chart belongs to a single datasource,
        and will have `params` that specifies which columns from that datasource are used.
        For example, if you need data from both 'cells' and 'genes' datasources, make sure to indicate this.
        
        Return a structured output describing a view configuration containing:
        - dataSources: A dictionary of datasource names to their configuration
        - initialCharts: A dictionary of datasource names to their chart configurations
        - title: The title of the view
        - param: A list of field specifications defining the data columns used by this chart
        """
        
        messages = [HumanMessage(content=analysis_prompt)]
        result = structured_llm.invoke(messages)
        return result  # type: ignore

    def generate_chart_structured(self, llm: ChatOpenAI, analysis: DataAnalysisResult, 
                                question: str, project: MDVProject) -> ChartGenerationResult:
        """Generate chart configuration and code using structured output"""
        
        # Create structured LLM for chart generation
        structured_llm = llm.with_structured_output(ChartGenerationResult)
        
        # Determine the best chart type
        chart_type = analysis.chart_types[0] if analysis.chart_types else "scatter_plot"
        
        # Parse selected columns to determine which datasources are involved
        datasources_involved = set()
        for column in analysis.selected_columns:
            if ':' in column:
                datasource_name = column.split(':')[0]
                datasources_involved.add(datasource_name)
            else:
                # If no datasource prefix, assume it's from the primary datasource
                datasources_involved.add(self.ds_name)
        
        # If no datasources were identified, use the primary one
        if not datasources_involved:
            datasources_involved = {self.ds_name}
        
        # Create chart configuration
        chart_config = create_chart_config(
            chart_type=chart_type,
            title=f"Visualization: {question}",
            legend=f"Chart showing {', '.join(analysis.selected_columns)}",
            params=analysis.selected_columns
        )
        
        # Generate Python code using RAG
        files_in_dir = os.listdir(project.dir)
        csv_file = None
        h5ad_file = None
        
        for file in files_in_dir:
            if file.endswith(".csv"):
                csv_file = file
            elif file.endswith(".h5ad"):
                h5ad_file = file
        
        if csv_file:
            path_to_data = os.path.join(project.dir, csv_file)
        elif h5ad_file:
            path_to_data = os.path.join(project.dir, h5ad_file)
        else:
            raise FileNotFoundError("No CSV or H5AD file found in the directory.")
        
        # Use the primary datasource name for RAG
        datasource_name = self.ds_name
        
        # Use RAG to generate code
        prompt_RAG = get_createproject_prompt_RAG(
            project, path_to_data, datasource_name, 
            str(analysis.selected_columns), question
        )
        
        prompt_RAG_template = PromptTemplate(
            template=prompt_RAG,
            input_variables=["context", "question"]
        )
        
        qa_chain = RetrievalQA.from_llm(
            llm=llm,
            prompt=prompt_RAG_template,
            retriever=retriever,
            return_source_documents=True,
        )
        
        output_qa = qa_chain.invoke({"query": question})
        result = output_qa["result"]
        
        # Prepare the code
        final_code = prepare_code(
            result,
            None,  # No longer using DataFrames for analysis
            project,
            self.log,
            modify_existing_project=True,
            view_name=question,
        )
        
        # Parse view name
        view_name = parse_view_name(final_code)
        if view_name is None:
            view_name = f"view_{hash(question)}"
        
        # Create chart generation result
        chart_result = ChartGenerationResult(
            chart_config=chart_config,
            python_code=final_code,
            view_name=view_name,
            dependencies=["mdvtools", "pandas", "scanpy"]
        )
        
        return chart_result

    def ask_question(self, chat_request: ChatRequest) -> AskQuestionResult:
        """
        Ask a question using structured analysis instead of PythonAstREPLTool
        """
        id = chat_request["id"]
        room = chat_request["room"]
        conversation_id = chat_request["conversation_id"]
        question = chat_request["message"]
        handle_error = chat_request["handle_error"]
        
        # Create socket API for this request
        socket_api = ChatSocketAPI(self.project, id, room, conversation_id)
        log = socket_api.log
        log(f"Asking question (structured): {question}")

        if question == "test error":
            raise Exception("testing error response as requested")
        
        # Ensure we have a conversation_id
        if not conversation_id:
            conversation_id = f"default_{id}"
        
        # Get or create memory for this conversation
        memory = self.get_or_create_memory(conversation_id)
        
        # Create local LLM instances with local callbacks for logging
        langchain_logging_handler = LangchainLoggingHandler(socket_api.logger)
        
        with time_block("b6: Initialising LLM for structured analysis"):
            structured_llm = ChatOpenAI(
                temperature=0.1, 
                model="gpt-4.1",
                callbacks=[langchain_logging_handler]
            )
        
        progress = 0
        
        # Initialize variables to avoid unbound variable errors in exception handling
        data_analysis = None
        chart_generation = None
        
        try:
            with time_block("b10a: Structured analysis"):
                socket_api.update_chat_progress(
                    "Structured data analysis", id, progress, 0
                )
                log(f"Analyzing data for: '{question}'")
                if self.init_error:
                    socket_api.update_chat_progress(
                        f"Agent initialisation error: {str(self.error_message)}", id, 100, 0
                    )
                    raise Exception(f"{str(self.error_message)[:500]}")
                logger.info(f"Question asked by user (structured): {question}")
            
            with time_block("b10b: Data analysis"):
                socket_api.update_chat_progress(
                    "Analyzing data...", id, progress, 40
                )
                progress += 40
                
                # Analyze data using structured output
                data_analysis = self.analyze_data_structured(
                    structured_llm, question, self.project
                )
                chat_debug_logger.info(f"Data Analysis Result: {data_analysis}")
            
            with time_block("b11: Chart generation"):
                socket_api.update_chat_progress(
                    "Generating chart...", id, progress, 30
                )
                progress += 30
                
                # Generate chart using structured output
                chart_generation = self.generate_chart_structured(
                    structured_llm, data_analysis, question, self.project
                )
                chat_debug_logger.info(f"Chart Generation Result: {chart_generation}")

            with time_block("b14: Execute code"):
                socket_api.update_chat_progress(
                    "Executing code...", id, progress, 30
                )
                progress += 30
                
                ok, stdout, stderr = execute_code(
                    chart_generation.python_code, open_code=False, log=log
                )
                chat_debug_logger.info(f"Code Execution Result - OK: {ok}")
                chat_debug_logger.info(f"Code Execution STDOUT:\n{stdout}")
                chat_debug_logger.info(f"Code Execution STDERR:\n{stderr}")

                if not ok:
                    raise Exception(f"Code execution failed: \n{stderr}")
                
                # Create structured result
                structured_result = StructuredQueryResult(
                    question=question,
                    data_analysis=data_analysis,
                    chart_generation=chart_generation,
                    execution_result={"stdout": stdout, "stderr": stderr, "success": ok},
                    error=None,
                    success=ok
                )
                
                # Set the view in the project with the generated chart
                # Determine which datasources are involved in this visualization
                datasources_involved = set()
                for column in data_analysis.selected_columns:
                    if ':' in column:
                        datasource_name = column.split(':')[0]
                        datasources_involved.add(datasource_name)
                    else:
                        # If no datasource prefix, assume it's from the primary datasource
                        datasources_involved.add(self.ds_name)
                
                # If no datasources were identified, use the primary one
                if not datasources_involved:
                    datasources_involved = {self.ds_name}
                
                # Create view configuration with appropriate panels
                view_config = {
                    "dataSources": {},
                    "initialCharts": {}
                }
                
                # Calculate panel widths based on number of datasources
                panel_width = 100 // len(datasources_involved)
                
                for i, ds_name in enumerate(datasources_involved):
                    view_config["dataSources"][ds_name] = {
                        "layout": "gridstack",
                        "panelWidth": panel_width
                    }
                    
                    # Add the chart to the appropriate datasource panel
                    # For now, put the chart in the first datasource panel
                    if i == 0:
                        view_config["initialCharts"][ds_name] = [chart_generation.chart_config.dict()]
                    else:
                        view_config["initialCharts"][ds_name] = []
                
                self.project.set_view(chart_generation.view_name, view_config)
                
            with time_block("b16: Log chat item"):
                final_code_updated = f"I ran some code for you (structured analysis):\n\n```python\n{chart_generation.python_code}\n```"
                log_chat_item(
                    self.project, question, 
                    {"structured_analysis": data_analysis.dict(), "chart_generation": chart_generation.dict()},
                    "Structured analysis approach", final_code_updated, conversation_id, chart_generation.view_name
                )
                log(final_code_updated)
                socket_api.update_chat_progress(
                    "Finished processing query (structured)", id, 100, 0
                )
                
                return {
                    "code": final_code_updated, 
                    "view_name": chart_generation.view_name, 
                    "error": False, 
                    "message": "Success (structured)"
                }
                
        except Exception as e:
            error_message = f"{str(e)[:500]}\n\n{traceback.format_exc()}"
            print(f"{error_message}")
            socket_api.update_chat_progress(
                f"Error: {error_message}", id, 100, 0
            )
            handle_error(e)
            
            # Create structured result with error information
            error_structured_result = StructuredQueryResult(
                question=question,
                data_analysis=data_analysis if data_analysis is not None else DataAnalysisResult(
                    selected_columns=[],
                    chart_types=[],
                    reasoning="Analysis failed due to error",
                    data_summary={}
                ),
                chart_generation=chart_generation if chart_generation is not None else ChartGenerationResult(
                    chart_config=create_chart_config("error", "Error", "Chart generation failed", []),
                    python_code="",
                    view_name="error",
                    dependencies=[]
                ),
                execution_result=None,
                error=error_message,
                success=False
            )
            
            return {
                "code": None, 
                "view_name": None, 
                "error": True, 
                "message": f"ERROR: {error_message}"
            } 