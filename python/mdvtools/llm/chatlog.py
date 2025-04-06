import datetime
import os
import json
from pathlib import Path
import logging
from typing import Any, Dict, List

from langchain_core.callbacks import BaseCallbackHandler
from langchain_core.messages import BaseMessage
from langchain_core.outputs import LLMResult


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
    print(f"Error checking log file exists: {e}")

def log_chat(output: Any, prompt_template: str, response: str):
    """
    output: result of invoke 'from langchain.chains import RetrievalQA'
    """
    print("Logging to json file...")
    context_information = output['source_documents']#output['context']
    context_information_metadata = [context_information[i].metadata for i in range(len(context_information))]
    context_information_metadata_url = [context_information_metadata[i]['url'] for i in range(len(context_information_metadata))]

    if file is not None:
        try:
            log_to_json(context_information_metadata_url, output['query'], prompt_template, response)#output['input'], prompt_template, response)
        except Exception as e:
            print(f"Error logging to json file: {e}")
    else:
        print("Json file does not exist. Skipping logging...")

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
