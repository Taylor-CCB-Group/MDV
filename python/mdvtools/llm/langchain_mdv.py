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
from .markdown_utils import create_suggested_questions_prompt

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
from .code_execution import execute_code
from .chatlog import LangchainLoggingHandler
from .mcp_tools import get_mcp_tools, init_mcp_tools, is_mcp_available

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
        
        # Initialize MCP tools with project reference
        try:
            init_mcp_tools(project)
            self.mcp_available = is_mcp_available()
            if self.mcp_available:
                logger.info("MCP-Bio server is available")
            else:
                logger.info("MCP-Bio server not available - MCP tools disabled")
        except Exception as e:
            logger.warning(f"Could not initialize MCP tools: {e}")
            self.mcp_available = False

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
            # Prepare dataframes for the agent
            if len(self.project.datasources) >= 2:
                self.df_list = [self.project.get_datasource_as_dataframe(ds['name']) for ds in self.project.datasources[:2]]
            else:
                self.df_list = [self.project.get_datasource_as_dataframe(self.project.datasources[0]['name'])]

            try:
                with time_block("b_suggest: Generating suggested questions"):
                    llm = ChatOpenAI(temperature=0.1, model="gpt-4o")
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
        Also includes MCP-Bio tools for bioinformatics analyses if available.
        """
        # Step 1: Create the Python REPL Tool
        python_tool = PythonAstREPLTool()
        
        if python_tool.globals is None:
            python_tool.globals = {}  # Ensure it's a dictionary
        
        # Make DataFrames available inside the REPL tool
        python_tool.globals.update(dfs) 
        
        # New fixes:
        python_tool.globals["list_globals"] = lambda: list(python_tool.globals.keys()) # type: ignore

        # Step 2: Gather all tools (Python REPL + MCP tools)
        all_tools = [python_tool]
        
        # Add MCP-Bio tools info
        if getattr(self, 'mcp_available', False):
            mcp_tools = get_mcp_tools()
            all_tools.extend(mcp_tools)
            # #region agent log
            log_debug("H1", "Tools registered with agent", {"tools": [t.name if hasattr(t, 'name') else str(t) for t in all_tools]})
            # #endregion
            mcp_tools_info = """

**Analysis Tools** (from MCP-Bio server):
- `list_analysis_tools`: See available bioinformatics analyses (UMAP, clustering, etc.)
- `run_umap`: Run UMAP dimensionality reduction
- `run_clustering`: Run Leiden clustering to group cells

For UMAP, clustering, or other bioinformatics analyses, use the analysis tools.
After running an analysis, tell the user to reload the page to see the new columns.
"""
            logger.info(f"Added {len(mcp_tools)} MCP tools to agent")
        else:
            mcp_tools_info = """

**Analysis Tools Status**: The MCP-Bio server is currently **OFFLINE** or unreachable at http://localhost:8000. 
If the user asks about available analysis tools (like UMAP, clustering, etc.), you must inform them that the MCP-Bio server is down and these tools are currently unavailable. 
**CRITICAL**: In this case, do NOT provide "fields" and "charts" lines in your output. Just explain that the server is offline. Do NOT try to generate data exploration code or charts as a substitute.
"""
            logger.info("MCP-Bio server not available - MCP tools disabled in prompt")

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

        # Step 4: Define the Agent Prompt (updated to include MCP tools)
        prompt_data_template = f"""You are an AI assistant for biological data analysis.

**Data Exploration Tools:**
You have access to the following Pandas DataFrames: {', '.join(dfs.keys())}. 
These are preloaded, so do not redefine them.

If the user is asking for a list of available tools or to run a specific analysis (like UMAP or clustering), use the Analysis Tools immediately. 
Otherwise, for data exploration questions, you must first run `df1.columns` and `df2.columns` to inspect available fields before answering. 
Use these to correct the column names mentioned by the user.
You must always invoke the PythonAstREPLTool to check the DataFrames columns and explore the values of the DataFrames.
Use `df.info()` or `df.index()`.
Before running any code, check available variables using `list_globals()`.
{mcp_tools_info}
""" + prompt_data

        prompt = ChatPromptTemplate.from_messages([
            ("system", prompt_data_template),
            MessagesPlaceholder(variable_name="chat_history"),
            ("human", prompt_data_template + "{input}"),
            ("ai", "{agent_scratchpad}"),
        ])

        # Step 5: Create the Agent with all tools
        agent = create_openai_functions_agent(llm, all_tools, prompt)

        # Step 6: Wrap in an Agent Executor (Finalized Agent)
        agent_executor = AgentExecutor(agent=agent, tools=all_tools, memory=memory, verbose=verbose, return_intermediate_steps=True)

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
        # #region agent log
        log_debug("H2", "ask_question called", {"question": chat_request.get("message")})
        # #endregion
        log = socket_api.log
        log(f"Asking question: {question}")

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
                model="gpt-4o",
                callbacks=[langchain_logging_handler]
            )
        
        with time_block("b7: Initialising LLM for agent"):
            dataframe_llm = ChatOpenAI(
                temperature=0.1, 
                model="gpt-4o",
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
        if mock_agent:
            ok, strdout, stderr = execute_code(
                'import mdvtools\nprint("mdvtools import ok")'
            )
            return {"code": f"mdvtools import ok? {ok}\n\nstdout: {strdout}\n\nstderr: {stderr}", "view_name": None, "error": False, "message": "Success"}
        
        with time_block("b10a: Agent invoking"):  # ~0.005% of time
            socket_api.update_chat_progress(
                "Agent invoking", id, progress, 0
            )
            log(f"Asking the LLM: '{question}'")
            if self.init_error:
                socket_api.update_chat_progress(
                    f"Agent initialisation error: {str(self.error_message)}", id, 100, 0
                )
                raise Exception(f"{str(self.error_message)[:500]}")
            logger.info(f"Question asked by user: {question}")
        try:
            with time_block("b10b: Pandas agent invoking"):  # ~31.4% of time
                socket_api.update_chat_progress(
                    "Pandas agent invoking...", id, progress, 31
                )
                progress += 31
                response = agent(question)
                # #region agent log
                steps = []
                if "intermediate_steps" in response:
                    for action, obs in response["intermediate_steps"]:
                        steps.append({"tool": action.tool, "input": action.tool_input, "output": str(obs)[:100]})
                log_debug("H4", "Agent execution completed", {
                    "output": response.get("output"), 
                    "steps": steps
                })
                # #endregion
                chat_debug_logger.info(f"Agent Response - output: {response['output']}")
            
            # If the response is from an MCP tool (using markers or checking intermediate steps)
            agent_output = response.get("output", "")
            mcp_tool_used = False
            if "intermediate_steps" in response:
                for action, _ in response["intermediate_steps"]:
                    # If any tool other than the python REPL was used, consider it an MCP/Analysis tool
                    if action.tool != "python_repl_ast":
                        mcp_tool_used = True
                        break
            
            # Stricter early exit for MCP-Bio related queries to avoid hallucinated charts
            is_mcp_query = "mcp-bio" in question.lower() or "analysis tool" in question.lower()
            if mcp_tool_used or (is_mcp_query and "mcp-bio" in agent_output.lower()) or any(marker in agent_output for marker in ["**Available Analysis Tools:**", "**UMAP Complete!**", "**Clustering Complete!**", "**UMAP Failed:**", "**Clustering Failed:**", "**Error**:", "MCP-Bio server is currently offline"]):
                # If the output is a JSON string (typical for read_h5ad), wrap it in a code block
                if agent_output.strip().startswith("{") and agent_output.strip().endswith("}"):
                    agent_output = f"```json\n{agent_output}\n```"
                
                # During debug, provide more info in chat
                steps_info = ""
                if "steps" in locals() and steps:
                    steps_info = "\n\n**Debug - Tool Steps:**\n"
                    for s in steps:
                        steps_info += f"- Tool: `{s['tool']}`\n  - Input: `{s['input']}`\n  - Output: `{s['output']}`\n"
                
                # Check for path-related errors and suggest troubleshooting
                if "File not found" in agent_output:
                    steps_info += "\n\n**Troubleshooting:**\n- The MCP server is on your host, while the app is in Docker.\n- We tried relative path (from mdv root) and absolute path (from container root).\n- Check if the MCP server was started in the same directory that contains the `2` folder."
                
                final_response = f"{agent_output}{steps_info}"
                log(final_response)
                socket_api.update_chat_progress("Finished processing analysis query", id, 100, 0)
                return {"code": final_response, "view_name": None, "error": False, "message": "Success"}

            with time_block("b11: RAG prompt preparation"):  # ~0.003% of time
                socket_api.update_chat_progress(
                    "RAG prompt preparation...", id, progress, 1
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
                    "Invoking RAG chain...", id, progress, 60
                )
                progress += 60

                match = re.search(r'charts\s+(.*)', response['output'])
                charts_part = match.group(1) if match else response['output']

                qa_chain = RetrievalQA.from_llm(
                    llm=code_llm, 
                    prompt=prompt_RAG_template, 
                    retriever=retriever,
                    return_source_documents=True,
                )
                output_qa = qa_chain.invoke({"query": charts_part})
                # Log the complete prompt after context injection
                if "source_documents" in output_qa:
                    retrieved_context = "\n".join(doc.page_content for doc in output_qa["source_documents"])
                    complete_prompt = prompt_RAG_template.format(context=retrieved_context, question=charts_part)
                    chat_debug_logger.info(f"=== Complete RAG Prompt (with injected context) ===\n{complete_prompt}\n=== End Complete Prompt ===")
                result = output_qa["result"]

            with time_block("b13: Prepare code"):  # <0.1% of time
                final_code = prepare_code(
                    result,
                    self.df,
                    self.project,
                    log,
                    modify_existing_project=True,
                    view_name=question,
                )

            chat_debug_logger.info(f"Prepared Code for Execution:\n{final_code}")
            chat_debug_logger.info(f"RAG output:\n{output_qa}")
            with time_block("b14: Execute code"):  # ~9% of time
                socket_api.update_chat_progress(
                    "Executing code...", id, progress, 9
                )
                progress += 9
                ok, stdout, stderr = execute_code(
                    final_code, open_code=False, log=log
                )
                chat_debug_logger.info(f"Code Execution Result - OK: {ok}")
                chat_debug_logger.info(f"Code Execution STDOUT:\n{stdout}")
                chat_debug_logger.info(f"Code Execution STDERR:\n{stderr}")

            if not ok:
                # Log code execution error
                raise Exception(f"Code execution failed: \n{stderr}")
            
            with time_block("b15: Parse view name"):
                # Parse view name from the code
                view_name = parse_view_name(final_code)
                if view_name is None:
                    raise Exception("Parsing view name failed")
                log(f"view_name: {view_name}")
                
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
                    f"{html.escape(explanation)}"
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
                    view_name=view_name
                )
                log(final_code_updated)
                socket_api.update_chat_progress(
                    "Finished processing query", id, 100, 0
                )
                # we want to know the view_name to navigate to as well... for now we do that in the calling code
                return {"code": final_code_updated, "view_name": view_name, "error": False, "message": "Success"}
        except Exception as e:
            # Log general error
            error_message = f"{str(e)[:500]}\n\n{traceback.format_exc()}"
            print(f"{error_message}")
            socket_api.update_chat_progress(
                f"Error: {error_message}", id, 100, 0
            )
            handle_error(e)
            return {"code": None, "view_name": None, "error": True, "message": f"ERROR: {error_message}"}

# #region agent log
import json
import time
def log_debug(hypothesis_id, message, data=None):
    log_entry = {
        "id": f"log_{int(time.time())}_{hypothesis_id}",
        "timestamp": int(time.time() * 1000),
        "location": "langchain_mdv.py",
        "message": message,
        "data": data or {},
        "sessionId": "debug-session",
        "runId": "run1",
        "hypothesisId": hypothesis_id
    }
    try:
        with open("/app/.cursor/debug.log", "a") as f:
            f.write(json.dumps(log_entry) + "\n")
    except Exception:
        pass
# #endregion
