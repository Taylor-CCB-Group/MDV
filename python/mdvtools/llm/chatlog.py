# https://stackoverflow.com/questions/39740632/python-type-hinting-without-cyclic-imports
from __future__ import annotations
import datetime
import os
import json
from pathlib import Path
import logging
from mdvtools.logging_config import get_logger

from typing import List, Optional, Dict, Any, TYPE_CHECKING

from langchain_core.callbacks import BaseCallbackHandler
from langchain_core.messages import BaseMessage
from langchain_core.outputs import LLMResult

from dataclasses import dataclass
# from mdvtools.websocket import log
from flask_socketio import SocketIO

if TYPE_CHECKING:
    from mdvtools.mdvproject import MDVProject


# all getting a bit messy and confusing in here - will likely be further refactored
logger = get_logger(__name__)
log = logger.info

@dataclass
class ChatLogItem:
    """Represents a chat log entry for a request and response"""
    context: str
    query: str
    prompt_template: str
    response: str
    timestamp: str
    conversation_id: Optional[str] = None
    view_name: Optional[str] = None
    error: Optional[bool] = None
    

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return {
            "context": self.context,
            "query": self.query,
            "prompt_template": self.prompt_template,
            "response": self.response,
            "timestamp": self.timestamp,
            "conversation_id": self.conversation_id,
            "view_name": self.view_name,
            "error": self.error,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ChatLogItem':
        """Create from dictionary (e.g. from JSON)"""
        return cls(
            context=data["context"],
            query=data["query"],
            prompt_template=data["prompt_template"],
            response=data["response"],
            timestamp=data["timestamp"],
            conversation_id=data.get("conversation_id"),
            view_name=data.get("view_name"),
            error=data.get("error"),
        )

class ChatLogger:
    """Handles chat logging functionality"""
    def __init__(self, log_file_path: str):
        self.log_file_path = Path(log_file_path)
        # note - not calling _ensure_log_file_exists here, 
        # it will currently be called when logging or reading logs
        # this is not necessarily the best design, but should avoid some current issues
        # and also I don't think we should be so eager to create the log file anyway.

    def _ensure_log_file_exists(self):
        """Ensure the log file exists and is properly initialized"""
        if not self.log_file_path.exists():
            self.log_file_path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.log_file_path, 'w') as f:
                json.dump([], f)
    def log_chat(self, item: ChatLogItem):
        """Log a chat item to the JSON file"""
        try:
            self._ensure_log_file_exists()
            # Read existing logs
            with open(self.log_file_path, 'r') as f:
                logs = json.load(f)
            
            # Append new log
            logs.append(item.to_dict())
            
            # Write back
            with open(self.log_file_path, 'w') as f:
                json.dump(logs, f, indent=4)
        except Exception as e:
            logging.error(f"Error logging chat: {e}")

    def get_logs(self) -> List[ChatLogItem]:
        """Get all chat logs"""
        try:
            self._ensure_log_file_exists()
            with open(self.log_file_path, 'r') as f:
                logs = json.load(f)
            return [ChatLogItem.from_dict(log) for log in logs]
        except Exception as e:
            logging.error(f"Error reading chat logs: {e}")
            return []

    def get_conversation_logs(self, conversation_id: str) -> List[ChatLogItem]:
        """Get logs for a specific conversation"""
        return [log for log in self.get_logs() if log.conversation_id == conversation_id] 

class ChatSocketAPI:
    """
    An instance of this class is created for each chat request.
    It will instantiate a Logger instance & SocketIOHandler which should be GCed when the request is finished.
    """
    def __init__(self, project: MDVProject, id: str, room: str, conversation_id: str):
        self.project = project
        from mdvtools.websocket import socketio
        if socketio is None:
            raise ValueError("SocketIO is not initialized")

        self.socketio = socketio
        # todo refactor event/room/namespace names 
        log_name = "chat"
        self.progress_name = "chat_progress"
        # we need to reevaluate so that 'logging' isn't sent to all clients.
        logger = logging.Logger(f"{log_name}_{project.id}_{id}")
        logger.propagate = False # avoid unintentional memory retention
        logger.setLevel(logging.INFO)
        self.project_namespace = f"/project/{project.id}"
        self.room = room
        handler = ChatSocketIOHandler(socketio, log_name, self.project_namespace, id, room)
        # todo - move more of this to another file, also add a handler to output to a file in project directory
        # probably don't need to have separate class for ChatSocketIOHandler, we could make this inherit logging.StreamHandler
        # would that be better, worse, or indifferent?
        logger.addHandler(handler)
        log_dir = os.path.join(project.dir, "logs")
        os.makedirs(log_dir, exist_ok=True)
        file_handler = logging.FileHandler(os.path.join(project.dir, f"logs/chat__{conversation_id}.log"))
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        file_handler.setLevel(logging.INFO)
        logger.addHandler(file_handler)
        self.logger = logger
        # self.log(f"ChatSocketAPI initialized for request {id} in room {room}")
        # self.update_chat_progress("ChatSocketAPI initialized", id, 0, 0)
        # from time import sleep
        # sleep(1)

    def log(self, msg: str):
        self.logger.info(msg)
    
    def update_chat_progress(self, message: str, id: str, progress: int, delta: int):
        """
        Send a message to the chat log and also update the progress bar.
        
        Args:
            message (str): the message to send to the chat log
            id (str): the id of the associated chat request
            progress (int): the progress value (%) to update the progress bar with
            delta (int): the expected cost of the current operation (%)
        """
        # we should descriminate which user to send this to... 
        # which implies that this instance should be associated...
        # or that we can map the chat request ID to a user ID.

        # I think simplest way to do this will be to use request.sid,
        # which implies that the `/chat` endpoint should be socket.io rather than REST.

        # to = chat_sid_map.get(id, None)
        # if to is None:
        #     log(f"Chat progress update for {id} but no associated user found, skipping.")
        #     return
        self.socketio.emit(self.progress_name, {
            "message": message, "id": id, "progress": progress, "delta": delta
        }, namespace=self.project_namespace, to=self.room)

class ChatSocketIOHandler(logging.StreamHandler):
    def __init__(self, socketio: SocketIO, event_name: str, namespace: str, id: str, room: str):
        super().__init__()
        log(f"handler initialized for event: {event_name}")
        self.socketio = socketio
        # todo - event_name vs namespace vs room refactor
        self.event_name = event_name
        self.namespace = namespace
        self.id = id
        self.room = room
        # def my_function_handler(data):
        #     log(f"{namespace}/{event_name}: {data}")
        # we could handle a cancel event handler here?
        # socketio.on_event(event_name, my_function_handler, namespace=namespace)

    def emit(self, record):
        """
        Emit a record - send it via socketio & also print it to the console.
        subject to change.
        """
        try:
            msg = self.format(record)
            log(f"[ {self.event_name} #{self.id} ] {msg}")
            #!!! to=self.id?
            self.socketio.emit(self.event_name, msg, namespace=self.namespace, to=self.room)
        except Exception:
            self.handleError(record)


mypath = os.path.dirname(__file__)
# this needs to be reviewed - causing trouble when running in different contexts
json_keyfile_path = os.path.join(mypath, "../../../chatlog.json")

# Function to ensure the JSON log file exists
def initialize_json_log():
    if not os.path.exists(json_keyfile_path):
        with open(json_keyfile_path, 'w') as file:
            json.dump([], file)  # Initialize with an empty list

# Function to log data to the JSON file
def log_to_json(context, prompt, prompt_template, response):
    # Ensure the log file exists
    initialize_json_log()

    # Prepare log entry
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
    log_entry = {
        "timestamp": timestamp,
        "context": context,
        "prompt": prompt,
        "prompt_template": prompt_template,
        "response": response
    }

    # Read the existing logs
    with open(json_keyfile_path, 'r') as file:
        logs = json.load(file)
    
    # Append the new log entry
    logs.append(log_entry)

    # Write back the updated logs to the file
    with open(json_keyfile_path, 'w') as file:
        json.dump(logs, file, indent=4)

try:
    file = Path(json_keyfile_path).exists()
except Exception as e:
    logging.error(f"Error checking log file exists: {e}")

class LangchainLoggingHandler(BaseCallbackHandler):
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.log = logger.info

    def on_chat_model_start(
        self, serialized: Dict[str, Any], messages: List[List[BaseMessage]], **kwargs
    ) -> None:
        self.log("Chat model started")

    def on_llm_end(self, response: LLMResult, **kwargs) -> None:
        self.log(f"Chat model ended, response: {response}")

    def on_chain_start(
        self, serialized: Dict[str, Any], inputs: Dict[str, Any], **kwargs
    ) -> None:
        self.log(f"Chain {serialized.get('name')} started")

    def on_chain_end(self, outputs: Dict[str, Any], **kwargs) -> None:
        self.log(f"Chain ended, outputs: {outputs}")

def log_chat_item(project, question, output, prompt_template, response, conversation_id, context: str | None, view_name: str | None, error: bool = False):
    """
    Log a chat interaction to the chat log file.
    Args:
        project: The MDVProject instance
        output: Result of invoke 'from langchain.chains import RetrievalQA' (can be None for errors)
        prompt_template: The template used for the prompt (can be empty for errors)
        response: The response generated (error message if error)
        conversation_id: ID to group messages from the same conversation
        context: Context of the response code generated which contains file names
        error: Whether this log is for an error
    """
    # Create a ChatLogger instance for this project
    chat_file = os.path.join(project.dir, "chat_log.json")
    chat_logger = ChatLogger(chat_file)
    if error or output is None:
        context = "[]"
        view_name = None
        prompt_template = prompt_template or ""
        

    # Ensure context is always a non-optional string for ChatLogItem
    context_str: str = context if context is not None else "[]"

    chat_item = ChatLogItem(
        context=context_str,
        query=question,
        prompt_template=prompt_template,
        response=response,
        timestamp=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        conversation_id=conversation_id,
        view_name=view_name,
        error=error
    )
    chat_logger.log_chat(chat_item)
