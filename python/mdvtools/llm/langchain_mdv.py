import time
import logging
import re
from contextlib import contextmanager
import os
import traceback
import html
from typing import List

# Code Generation using Retrieval Augmented Generation + LangChain
from mdvtools.llm.chat_protocol import AskQuestionResult, ChatRequest, ProjectChatProtocol
from mdvtools.mdvproject import MDVProject

from langchain_openai import ChatOpenAI
from langchain_openai import OpenAIEmbeddings
from langchain.schema.document import Document
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import FAISS
from langchain.text_splitter import Language
from langchain_core.pydantic_v1 import BaseModel, Field
from langchain.schema import HumanMessage
from mdvtools.markdown_utils import create_suggested_questions_prompt

# from langchain.prompts import PromptTemplate

from dotenv import load_dotenv
from mdvtools.llm.chatlog import ChatSocketAPI, LangchainLoggingHandler, log_chat_item

# packages for custom langchain agent
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain_experimental.tools.python.tool import PythonAstREPLTool
from langchain.agents import create_openai_functions_agent, AgentExecutor
from langchain.chains import LLMChain
from .local_files_utils import crawl_local_repo, extract_python_code_from_py, extract_python_code_from_ipynb
from .templates import get_createproject_prompt_RAG, prompt_data
from .code_manipulation import parse_view_name, prepare_code, extract_explanation_from_response
from .column_field_resolve import normalize_view_chart_params, prune_view_charts_with_invalid_params
from .verification import build_verification_summary
from .datasource_roles import collect_wrapper_subgroup_keys_for_project, infer_datasource_roles
from .code_execution import execute_code
from .chat_preview import format_stdout_for_chat
from .chatlog import LangchainLoggingHandler
from .chat_client_refresh import client_needs_refresh_after_chat
from .execution_progress import (
    attach_failed_source_context,
    ProgressEvent,
    ProgressThrottler,
    build_heartbeat_event,
    friendly_subprocess_failure_message,
    infer_progress_event_from_output,
    parse_explicit_progress_line,
)
from .preflight_flow import preflight_with_single_retry

# packages for memory
from langchain.chains import create_history_aware_retriever, create_retrieval_chain
from langchain.memory import ConversationBufferMemory

from langchain.prompts import PromptTemplate
from langchain.chains import RetrievalQA

import matplotlib

# create logger
logger = logging.getLogger('timing_results')
logger.setLevel(logging.DEBUG)

# create file handler and set level to INFO
file_handler = logging.FileHandler('timing_results.log')
file_handler.setLevel(logging.INFO)
logger.addHandler(file_handler)
logger.info(time.asctime())


# Create new file handler
# create a separate logger for chat debugging
chat_debug_logger = logging.getLogger('chat_debug')
chat_debug_logger.setLevel(logging.DEBUG)

# create file handler for chat_debug
chat_debug_handler = logging.FileHandler('chat_debug.log')
chat_debug_handler.setLevel(logging.DEBUG)

# optional: formatter (re-use if you already have one)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
chat_debug_handler.setFormatter(formatter)

# attach the handler
chat_debug_logger.addHandler(chat_debug_handler)


_CHART_CLASS_RE = re.compile(
    r"\b(DotPlot|ScatterPlot|HeatmapPlot|HistogramPlot|BoxPlot|ViolinPlot|"
    r"DensityScatterPlot|ScatterPlot3D|RowChart|StackedRowChart|PieChart|RingChart|"
    r"AbundanceBoxPlot|MultiLinePlot|TablePlot|WordcloudPlot|SankeyPlot|"
    r"SelectionDialogPlot|RowSummaryBox|TextBox)\s*\("
)


def _is_text_table_intent(question: str) -> bool:
    q = question.lower()
    text_table_keywords = [
        "predict",
        "annotate",
        "annotation",
        "label",
        "list",
        "table",
        "summarize",
        "summary",
        "top ",
        "rank",
        "mapping",
        "which cell type",
        "cell type",
    ]
    chart_intent_keywords = [
        "plot",
        "chart",
        "visualize",
        "scatter",
        "heatmap",
        "umap",
        "tsne",
        "violin",
        "histogram",
        "dotplot",
    ]
    has_text_table = any(k in q for k in text_table_keywords)
    has_explicit_chart = any(k in q for k in chart_intent_keywords)
    return has_text_table and not has_explicit_chart


def _extract_chart_classes(code: str) -> list[str]:
    seen: list[str] = []
    for m in _CHART_CLASS_RE.finditer(code):
        c = m.group(1)
        if c not in seen:
            seen.append(c)
    return seen


def _is_text_table_only_code(code: str) -> bool:
    chart_classes = _extract_chart_classes(code)
    if not chart_classes:
        return True
    allowed = {"TextBox", "TablePlot", "SelectionDialogPlot", "RowSummaryBox"}
    return all(c in allowed for c in chart_classes)


def _build_rag_retrieval_query(user_question: str, charts_part: str) -> str:
    """Combine the user question with agent-suggested chart types for RAG retrieval."""
    question = (user_question or "").strip()
    charts = (charts_part or "").strip()
    if question and charts:
        return f"User question: {question}\nSuggested chart types: {charts}"
    return question or charts


class SuggestedQuestions(BaseModel):
    """A list of suggested questions to ask about the data."""
    questions: List[str] = Field(
        description="A list of 5 suggested questions based on the provided data context. Each question should be a single, complete sentence."
    )

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
# if it isn't set, we raise an error more explicitly to log it
if not os.getenv("OPENAI_API_KEY"):
    print("OPENAI_API_KEY is not set in the environment variables. Please set it if LLM integration is required.")
    raise ValueError("OPENAI_API_KEY is not set in the environment variables. Please set it if LLM integration is required.")

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


class ProjectChat(ProjectChatProtocol):
    def __init__(self, project: MDVProject):
        self.project = project
        # could come from a config file - was passed as argument previously
        # but I'm reducing how far this prototype reaches into wider code
        self.welcome = (
            "Hello, I'm an AI assistant that has access to the data in this project "
            "and is designed to help build views for visualising it. What can I help you with?"
        )

        # rather than assign socket_api.logger to self.log, we can distinguish this global logger
        # from the one used during a chat request.
        self.log = logger.info
        
        # Store conversation memories for persistence across requests
        self.conversation_memories = {}  # Dict to store ConversationBufferMemory for each conversation
        self.suggested_questions: list[str] = []
        if len(project.datasources) == 0:
            raise ValueError("The project does not have any datasources")
        elif len(project.datasources) > 1: # remove? or make it == 1 ?
            # tempting to think we should just give the agent all of the datasources here
            self.log("The project has more than one datasource, only the first one will be used")
            self.ds_name1 = project.datasources[1]['name'] # maybe comment this out?
            self.df1 = project.get_datasource_as_dataframe(self.ds_name1) # and maybe comment this out?

        self.ds_name = project.datasources[0]['name']
        try:
            # raise ValueError("test error")
            self.df = None
            # Prepare dataframes for the agent.
            # ChatMDV historically assumed df1=cells and df2=genes. Many projects instead use
            # wrapper-capable expression datasources (e.g. rna/protein) linked via rows-as-columns.
            roles = infer_datasource_roles(self.project)
            self.roles = roles
            df1 = self.project.get_datasource_as_dataframe(roles.obs_datasource)
            primary_expr = roles.preferred_expression()
            if primary_expr is not None:
                df2 = self.project.get_datasource_as_dataframe(primary_expr.datasource_name)
                self.df_list = [df1, df2]
            else:
                self.df_list = [df1]

            try:
                with time_block("b_suggest: Generating suggested questions"):
                    llm = ChatOpenAI(temperature=0.1, model="gpt-4.1")
                    structured_llm = llm.with_structured_output(SuggestedQuestions)
                    prompt_text = create_suggested_questions_prompt(self.project)
                    messages = [HumanMessage(content=prompt_text)]
                    result = structured_llm.invoke(messages)
                    self.suggested_questions = result.questions # type: ignore
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

    def create_custom_pandas_agent(self, llm, dfs: dict, prompt_data, memory, verbose=False):
        """
        Creates a LangChain agent that can interact with Pandas DataFrames using a Python REPL tool.
        """
        # Step 1: Create the Python REPL Tool
        python_tool = PythonAstREPLTool()
        
        if python_tool.globals is None:
            python_tool.globals = {}  # Ensure it's a dictionary
        
        # Make DataFrames available inside the REPL tool
        python_tool.globals.update(dfs) 
        
        # New fixes:
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
        
        # > LangChainDeprecationWarning: The class `LLMChain` was deprecated in LangChain 0.1.17 and will be removed in 1.0. 
        # Use RunnableSequence, e.g., `prompt | llm` instead.
        contextualize_chain = LLMChain(llm=llm, prompt=contextualize_prompt, memory=memory)

        # Step 3: Define the Agent Prompt
        role_hint = ""
        try:
            roles = getattr(self, "roles", None)
            if roles is not None:
                exprs = roles.expressions or []
                if exprs:
                    expr_lines = "\n".join(
                        f"- df2 maps to expression datasource '{exprs[0].datasource_name}' "
                        f"(feature names in column '{exprs[0].name_column}', default subgroup '{exprs[0].subgroup_key}')"
                    )
                else:
                    expr_lines = "- No rows-as-columns expression datasource detected for this project."
                role_hint = (
                    "\n\nDatasource roles:\n"
                    f"- df1 maps to obs datasource '{roles.obs_datasource}'\n"
                    f"{expr_lines}\n"
                    "- For any chart `params` on the **feature table datasource** (df2), use **Field ID** strings from "
                    "`df2.columns` / project metadata for **that** datasource (e.g. `gene_ids` when listed for `genes`); "
                    "do not use those ids as `params` on charts bound to **cells** unless `cells` lists the same Field ID.\n"
                    "- Do not pair Scanpy `rank_genes_groups` tables with DotPlot/Heatmap on cells using wrapper `params` "
                    "as a substitute for that table (see ChatMDV \"Marker ranking vs DotPlot\" in RAG); optional "
                    "`add_datasource('chat_rank_genes_result', ...)` for an in-view table.\n"
                    "- Do not pair other Scanpy summary tables with MDV wrapper-based expression charts for the "
                    "same metric unless values are guaranteed identical; prefer one pipeline (see ChatMDV viz consistency "
                    "in RAG).\n"
                )
        except Exception:
            role_hint = ""

        prompt_data_template = f"""You have access to the following Pandas DataFrames: 
        {', '.join(dfs.keys())}. These are preloaded, so do not redefine them.
        {role_hint}
        Before answering any user question, you must first run `df1.columns` and `df2.columns` to inspect available fields. 
        Use these to correct the column names mentioned by the user.
        Before writing any `project.get_datasource_as_dataframe(datasource_name, columns=[...])` call, you must verify the datasource schema and only request columns that exist on that datasource.
        Never assume a `name` column exists on `cells`; if the task asks for genes/markers, use the expression datasource or explicit marker ranking outputs instead of `cells.name`.
        You must always invoke the PythonAstREPLTool to check the DataFrames columns and explore the values of the DataFrames.
        Use `df.info()` or `df.index()`.
        Before running any code, check available variables using `list_globals()`.""" + prompt_data

        prompt = ChatPromptTemplate.from_messages([
            ("system", prompt_data_template),
            MessagesPlaceholder(variable_name="chat_history"),
            ("human", prompt_data_template + "{input}"),
            ("ai", "{agent_scratchpad}"),
        ])

        # Step 4: Create the Agent
        agent = create_openai_functions_agent(llm, [python_tool], prompt)

        # Step 5: Wrap in an Agent Executor (Finalized Agent)
        agent_executor = AgentExecutor(agent=agent, tools=[python_tool], memory=memory, verbose=verbose, return_intermediate_steps=True)

        # Step 6: Wrapper Function to Use Contextualization and Preserve Memory
        def agent_with_contextualization(question):
            standalone_question = contextualize_chain.run(input=question)
            # Point 1: Log reformulation
            chat_debug_logger.info(f"Original Question: {question}")
            chat_debug_logger.info(f"Standalone Reformulated Question: {standalone_question}")
            # Point 2: Log what you're sending to the agent
            chat_debug_logger.info(f"Sending to agent_executor with input: {standalone_question}")
            response = agent_executor.invoke({"input": standalone_question})
            # Point 3: Log agent output
            chat_debug_logger.info(f"Agent Raw Response: {response}")
            memory.save_context({"input": question}, {"output": response.get("output", str(response))})
            return response

        return agent_with_contextualization

    def ask_question(self, chat_request: ChatRequest) -> AskQuestionResult:
        """
        Ask a question, generate code to answer it, execute the code...
        If the question is "test error", we raise an error to test the error handling.
        """
        id = chat_request["id"]
        room = chat_request["room"]
        conversation_id = chat_request["conversation_id"]
        question = chat_request["message"]
        handle_error = chat_request["handle_error"]
        
        # Create socket API for this request
        socket_api = ChatSocketAPI(self.project, id, room, conversation_id)
        log = chat_debug_logger.info
        chat_debug_logger.info("Asking question: %s", question)

        if question == "test error":
            raise Exception("testing error response as requested")
        
        # Ensure we have a conversation_id
        if not conversation_id:
            conversation_id = f"default_{id}"
        
        # Get or create memory for this conversation
        memory = self.get_or_create_memory(conversation_id)
        
        # Create local LLM instances with local callbacks for logging
        langchain_logging_handler = LangchainLoggingHandler(socket_api.logger)
        
        with time_block("b6: Initialising LLM for RAG"):
            code_llm = ChatOpenAI(
                temperature=0.1, 
                model="gpt-4.1",
                callbacks=[langchain_logging_handler]
            )
        
        with time_block("b7: Initialising LLM for agent"):
            dataframe_llm = ChatOpenAI(
                temperature=0.1, 
                model="gpt-4.1",
                callbacks=[langchain_logging_handler]
            )
        # there is a risk that these are out of date by the time we get here...
        # but it's a bit of an expensive operation, so avoiding it for now (in most cases we shouldn't need to refer to them)
        df_list = self.df_list
        
        # Create agent with local LLM and memory
        with time_block("b8: Pandas agent creating"):
            agent = self.create_custom_pandas_agent(
                dataframe_llm, 
                {"df1": df_list[0], "df2": df_list[1] if len(df_list) > 1 else df_list[0]}, 
                prompt_data, 
                memory,
                verbose=True
            )
        
        progress = 0
        total_steps = 6

        def _step_message(step: int, text: str) -> str:
            return f"Step {step}/{total_steps}: {text}"
        if mock_agent:
            ok, strdout, stderr = execute_code(
                'import mdvtools\nprint("mdvtools import ok")'
            )
            return {
                "code": f"mdvtools import ok? {ok}\n\nstdout: {strdout}\n\nstderr: {stderr}",
                "view_name": None,
                "error": False,
                "message": "Success",
                "verification": None,
                "data_preview": format_stdout_for_chat(strdout),
                "needs_refresh": False,
            }
        
        with time_block("b10a: Agent invoking"):  # ~0.005% of time
            socket_api.update_chat_progress(
                _step_message(1, "Understanding your request"), id, progress, 0, step_index=1, step_total=total_steps
            )
            chat_debug_logger.info("Asking the LLM: %r", question)
            if self.init_error:
                socket_api.update_chat_progress(
                    f"Agent initialisation error: {str(self.error_message)}", id, 100, 0
                )
                raise Exception(f"{str(self.error_message)[:500]}")
            logger.info(f"Question asked by user: {question}")
        try:
            with time_block("b10b: Pandas agent invoking"):  # ~31.4% of time
                socket_api.update_chat_progress(
                    _step_message(2, "Inspecting datasources and available fields"),
                    id,
                    progress,
                    31,
                    step_index=2,
                    step_total=total_steps,
                )
                progress += 31
                socket_api.update_chat_progress(
                    _step_message(3, "Generating runnable analysis code (this can take up to 1-2 minutes)"),
                    id,
                    progress,
                    0,
                    step_index=3,
                    step_total=total_steps,
                )
                response = agent(question)
                chat_debug_logger.info(f"Agent Response - output: {response['output']}")
            
            with time_block("b11: RAG prompt preparation"):  # ~0.003% of time
                socket_api.update_chat_progress(
                    _step_message(3, "Preparing analysis context"),
                    id,
                    progress,
                    1,
                    step_index=3,
                    step_total=total_steps,
                )
                # List all files in the directory
                files_in_dir = os.listdir(self.project.dir)

                # Initialize variables
                csv_file = None
                h5ad_file = None

                # Identify the CSV or H5AD file (optional)
                # subject to review...
                for file in files_in_dir:
                    if file.endswith(".csv"):
                        csv_file = file
                    elif file.endswith(".h5ad"):
                        h5ad_file = file

                # Determine the path_to_data (optional)
                if h5ad_file:
                    path_to_data = os.path.join(self.project.dir, h5ad_file)
                elif csv_file:
                    path_to_data = os.path.join(self.project.dir, csv_file)
                else:
                    path_to_data = ""  # fallback: no external file; operate on existing project datasources

                datasource_names = [ds['name'] for ds in self.project.datasources[:2]]  # Get names of up to 2 datasources

                #!!!!!! for now, assuming there will be an anndata.h5ad file in the project directory and will fail ungacefully if there isn't!!!!
                # we pass a reference to the actual project object and let figuring out the path be an internal implementation detail...
                # this should be more robust, and also more flexible in terms of what reasoning this method may be able to do internally in the future
                # I appear to have an issue though - the configuration of the devcontainer doesn't flag whether or not the thing we're passing is the right type
                # and the assert in the function is being triggered even though it should be fine
                prompt_RAG = get_createproject_prompt_RAG(self.project, path_to_data, datasource_names[0], response['output'], response['input'])
                chat_debug_logger.info(f"=== RAG Base Prompt (before context injection) ===\n{prompt_RAG}\n=== End Base Prompt ===")

                prompt_RAG_template = PromptTemplate(
                     template=prompt_RAG,
                     input_variables=["context", "question"])

            with time_block("b12: RAG chain"):  # ~60% of time
                socket_api.update_chat_progress(
                    _step_message(3, "Generating runnable analysis code"),
                    id,
                    progress,
                    60,
                    step_index=3,
                    step_total=total_steps,
                )
                progress += 60

                match = re.search(r'charts\s+(.*)', response['output'])
                charts_part = match.group(1) if match else response['output']
                rag_retrieval_query = _build_rag_retrieval_query(question, charts_part)

                qa_chain = RetrievalQA.from_llm(
                    llm=code_llm, 
                    prompt=prompt_RAG_template, 
                    retriever=retriever,
                    return_source_documents=True,
                )
                output_qa = qa_chain.invoke({"query": rag_retrieval_query})
                # Log the complete prompt after context injection
                if "source_documents" in output_qa:
                    retrieved_context = "\n".join(doc.page_content for doc in output_qa["source_documents"])
                    complete_prompt = prompt_RAG_template.format(
                        context=retrieved_context, question=rag_retrieval_query
                    )
                    chat_debug_logger.info(f"=== Complete RAG Prompt (with injected context) ===\n{complete_prompt}\n=== End Complete Prompt ===")
                result = output_qa["result"]

            with time_block("b13: Prepare code"):  # <0.1% of time
                final_code = prepare_code(
                    result,
                    self.df,
                    self.project,
                    chat_debug_logger.info,
                    modify_existing_project=True,
                    view_name=question,
                )

                if _is_text_table_intent(question) and _is_text_table_only_code(final_code):
                    chat_debug_logger.info(
                        "Detected text/table-first intent with chart-only output; regenerating with chat-first instruction."
                    )
                    chat_first_query = (
                        f"{rag_retrieval_query}\n\n"
                        "IMPORTANT: This request is text/table-first. "
                        "Prefer chat explanation and concise bounded stdout (`print(...)`) instead of creating "
                        "TextBox/TablePlot/SelectionDialog-only outputs unless the user explicitly requested a chart."
                    )
                    output_qa_retry = qa_chain.invoke({"query": chat_first_query})
                    result_retry = output_qa_retry["result"]
                    final_code_retry = prepare_code(
                        result_retry,
                        self.df,
                        self.project,
                        chat_debug_logger.info,
                        modify_existing_project=True,
                        view_name=question,
                    )
                    if not _is_text_table_only_code(final_code_retry):
                        output_qa = output_qa_retry
                        result = result_retry
                        final_code = final_code_retry
                    else:
                        # Keep retry result anyway when it remains text/table-only;
                        # downstream chat explanation/stdout still carries the answer.
                        output_qa = output_qa_retry
                        result = result_retry
                        final_code = final_code_retry

            with time_block("b13b: Preflight validate and retry"):
                datasource_fields: dict[str, set[str]] = {}
                for ds in self.project.datasources or []:
                    if not isinstance(ds, dict):
                        continue
                    ds_name = ds.get("name")
                    cols = ds.get("columns")
                    if not isinstance(ds_name, str) or not isinstance(cols, list):
                        continue
                    datasource_fields[ds_name] = {
                        str(c.get("field"))
                        for c in cols
                        if isinstance(c, dict) and c.get("field")
                    }

                def _regenerate_for_preflight(issue_text: str) -> str:
                    retry_query = (
                        f"{rag_retrieval_query}\n\n"
                        "Fix this code preflight validation failure and return corrected runnable code only.\n"
                        "If issues mention set_x_axis or set_y_axis on BoxPlot, ViolinPlot, ScatterPlot, or DotPlot, "
                        "replace with set_axis_properties('x', {...}) and set_axis_properties('y', {...}) only.\n"
                        "On ScatterPlot use set_filter, not set_on_filter (set_on_filter is for ScatterPlot3D only).\n"
                        f"Preflight issues:\n{issue_text}"
                    )
                    output_qa_retry = qa_chain.invoke({"query": retry_query})
                    result_retry = output_qa_retry["result"]
                    return prepare_code(
                        result_retry,
                        self.df,
                        self.project,
                        chat_debug_logger.info,
                        modify_existing_project=True,
                        view_name=question,
                    )

                allowed_wrapper_subgroup_keys = collect_wrapper_subgroup_keys_for_project(
                    self.project
                )
                final_code, preflight_meta = preflight_with_single_retry(
                    initial_code=final_code,
                    regenerate_once=_regenerate_for_preflight,
                    log=chat_debug_logger.info,
                    datasource_fields=datasource_fields,
                    allowed_wrapper_subgroup_keys=allowed_wrapper_subgroup_keys,
                )
                chat_debug_logger.info("Preflight metadata: %s", preflight_meta)

            chat_debug_logger.info(f"Prepared Code for Execution:\n{final_code}")
            chat_debug_logger.info(f"RAG output:\n{output_qa}")
            with time_block("b14: Execute code"):  # ~9% of time
                socket_api.update_chat_progress(
                    _step_message(4, "Running analysis on your data"),
                    id,
                    progress,
                    9,
                    step_index=4,
                    step_total=total_steps,
                )
                progress += 9
                progress_throttler = ProgressThrottler(min_interval_seconds=3.0)
                current_stage: str | None = None

                def emit_progress_event(event: ProgressEvent) -> None:
                    nonlocal current_stage
                    if event.stage:
                        current_stage = event.stage
                    if event.progress == 0 and progress > 0:
                        event = ProgressEvent(
                            message=event.message,
                            progress=progress,
                            delta=event.delta,
                            stage=event.stage,
                            step_index=event.step_index,
                            step_total=event.step_total,
                            eta_seconds=event.eta_seconds,
                            elapsed_seconds=event.elapsed_seconds,
                            source=event.source,
                        )
                    if not progress_throttler.should_emit(event):
                        return
                    socket_api.update_chat_progress(
                        event.message,
                        id,
                        event.progress,
                        event.delta,
                        stage=event.stage,
                        step_index=event.step_index,
                        step_total=event.step_total,
                        eta_seconds=event.eta_seconds,
                        elapsed_seconds=event.elapsed_seconds,
                        source=event.source,
                    )

                def handle_output_line(_source: str, line: str) -> None:
                    explicit = parse_explicit_progress_line(line)
                    if explicit is not None:
                        emit_progress_event(explicit)
                        return
                    inferred = infer_progress_event_from_output(line)
                    if inferred is not None:
                        emit_progress_event(inferred)

                def handle_heartbeat(elapsed_seconds: float) -> None:
                    emit_progress_event(
                        build_heartbeat_event(
                            elapsed_seconds,
                            stage=current_stage,
                            step_index=4,
                            step_total=total_steps,
                        )
                    )

                ok, stdout, stderr = execute_code(
                    final_code,
                    open_code=False,
                    log=chat_debug_logger.info,
                    on_output_line=handle_output_line,
                    on_heartbeat=handle_heartbeat,
                )
                chat_debug_logger.info(f"Code Execution Result - OK: {ok}")
                chat_debug_logger.info(f"Code Execution STDOUT:\n{stdout}")
                chat_debug_logger.info(f"Code Execution STDERR:\n{stderr}")

            if not ok:
                # Preserve full diagnostics in backend logs, but surface friendly messaging for kill-style exits.
                chat_debug_logger.error("Code execution failed diagnostics:\n%s", stderr)
                friendly = friendly_subprocess_failure_message(stderr or "")
                if friendly is not None:
                    raise Exception(
                        attach_failed_source_context(
                            friendly,
                            stderr or "",
                            fallback_code=final_code,
                        )
                    )
                raise Exception(f"Code execution failed: \n{stderr}")

            data_preview_text = format_stdout_for_chat(stdout)

            with time_block("b15: Parse view name"):
                socket_api.update_chat_progress(
                    _step_message(5, "Validating and organizing results"),
                    id,
                    progress,
                    0,
                    step_index=5,
                    step_total=total_steps,
                )
                # Parse view name from the code
                view_name = parse_view_name(final_code)
                if view_name is None:
                    raise Exception("Parsing view name failed")
                log(f"view_name: {view_name}")

            with time_block("b15c: Normalize chart params (name to field id)"):
                try:
                    view_obj = self.project.get_view(view_name)
                    if view_obj:
                        ic = view_obj.get("initialCharts")
                        if isinstance(ic, dict):
                            valid_names = {
                                str(d["name"])
                                for d in (self.project.datasources or [])
                                if isinstance(d, dict) and d.get("name")
                            }
                            for ds_name in ic:
                                if ds_name not in valid_names:
                                    chat_debug_logger.warning(
                                        "View %s initialCharts references unknown datasource %r (have: %s)",
                                        view_name,
                                        ds_name,
                                        sorted(valid_names),
                                    )
                        normalized = normalize_view_chart_params(
                            view_obj, self.project
                        )
                        pruned, n_dropped = prune_view_charts_with_invalid_params(
                            normalized, self.project
                        )
                        if n_dropped:
                            log(
                                f"pruned {n_dropped} chart(s) with param tokens not in datasource metadata"
                            )
                        self.project.set_view(view_name, pruned)
                except Exception as norm_ex:
                    chat_debug_logger.warning(
                        "normalize_view_chart_params failed for view %s: %s",
                        view_name,
                        norm_ex,
                        exc_info=True,
                    )

            verification_text = ""
            with time_block("b15b: Verification summary"):
                try:
                    verification_text = build_verification_summary(
                        self.project, final_code, view_name
                    )
                except Exception as ver_ex:
                    chat_debug_logger.warning(
                        "build_verification_summary failed for view %s: %s",
                        view_name,
                        ver_ex,
                        exc_info=True,
                    )

            with time_block("b16: Log chat item"):
                 # Extract the explanation section from the LLM's response (removing code blocks)
                explanation = extract_explanation_from_response(output_qa["result"])
                # Extract the context information from the response
                context_information = output_qa['source_documents']
                context_information_metadata = [context_information[i].metadata for i in range(len(context_information))]
                context_information_metadata_url = [context_information_metadata[i]['url'] for i in range(len(context_information_metadata))]
                context_information_metadata_name = [s for s in context_information_metadata_url]
                context = str(context_information_metadata_name)
                context_files = (
                    "<br><br>"
                    "The context used to generate the above code has been augmented by the following files:\n\n"
                    # Prevent XSS by adding html.escape (Suggested by coderabbit)
                    + "\n".join(f"- `{html.escape(name)}`" for name in context_information_metadata_name)
                    + "\n\n"
                )
                final_code_updated = (
                    "I ran some code for you:\n\n"
                    "```python\n"
                    f"{final_code}\n"
                    "```\n\n"
                    "<br><br>"
                    f"{explanation}"
                    "\n\n"
                    f"{context_files}"
                )
                # Log successful code execution
                log_chat_item(
                    project=self.project,
                    question=question,
                    output=output_qa,
                    prompt_template=prompt_RAG,
                    response=final_code_updated,
                    conversation_id=conversation_id,
                    context=context,
                    view_name=view_name,
                    verification=verification_text or None,
                    data_preview=data_preview_text,
                )
                log(final_code_updated)
                socket_api.update_chat_progress(
                    _step_message(6, "Done — results are ready"),
                    id,
                    100,
                    0,
                    step_index=6,
                    step_total=total_steps,
                )
                # we want to know the view_name to navigate to as well... for now we do that in the calling code
                return {
                    "code": final_code_updated,
                    "view_name": view_name,
                    "error": False,
                    "message": "Success",
                    "verification": verification_text or None,
                    "data_preview": data_preview_text,
                    "needs_refresh": client_needs_refresh_after_chat(final_code),
                }
        except Exception as e:
            # Log general error
            error_message = f"{str(e)[:500]}\n\n{traceback.format_exc()}"
            print(f"{error_message}")
            socket_api.update_chat_progress(
                f"Analysis failed: {str(e)[:200]}",
                id,
                100,
                0,
                step_index=4,
                step_total=total_steps,
            )
            handle_error(e)
            return {
                "code": None,
                "view_name": None,
                "error": True,
                "message": f"ERROR: {error_message}",
                "verification": None,
                "data_preview": None,
                "needs_refresh": False,
            }
