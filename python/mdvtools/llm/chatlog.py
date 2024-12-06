import datetime
import os
import json
from typing import Any
from pathlib import Path


mypath = os.path.dirname(__file__)
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
    context_information = output['source_documents']
    context_information_metadata = [context_information[i].metadata for i in range(len(context_information))]
    context_information_metadata_url = [context_information_metadata[i]['url'] for i in range(len(context_information_metadata))]

    if file is not None:
        log_to_json(context_information_metadata_url, output['query'], prompt_template, response)
    else:
        print("Json file does not exist. Skipping logging...")
    