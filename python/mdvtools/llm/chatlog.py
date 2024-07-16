import datetime
# from google.oauth2 import service_account
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import os
from typing import Any

mypath = os.path.dirname(__file__)
json_keyfile_path = os.path.join(mypath, "key.json")
sheet_name = "Testing for RAG context"  # TODO: Update this with your sheet's name


# Function to initialize the Google Sheets connection
def init_google_sheet(json_keyfile_path, sheet_name):
    scopes = [
        "https://www.googleapis.com/auth/spreadsheets",
        "https://www.googleapis.com/auth/drive",
    ]
    credentials = ServiceAccountCredentials.from_json_keyfile_name(
        json_keyfile_path, scopes
    )
    gc = gspread.authorize(credentials)
    sheet = gc.open(sheet_name).sheet1  # Opens the first sheet in the spreadsheet
    return sheet


# Function to log data to the Google Sheet
def log_to_google_sheet(sheet: gspread.worksheet, context: str, prompt: str, prompt_template: str, response: str):
    timestamp = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d %H:%M:%S")
    data = [timestamp, context, prompt, prompt_template, response]
    sheet.append_row(data)  # Appends a row with the prompt and response

try:
    sheet = init_google_sheet(json_keyfile_path, sheet_name)
except Exception as e:
    print(f"Error initializing Google Sheet: {e}")


def log_chat(output: Any, prompt_template: str, response: str):
    """
    output: result of invoke 'from langchain.chains import RetrievalQA'
    """
    print("Logging to Google Sheet...")
    context_information = output['source_documents']
    context_information_metadata = [context_information[i].metadata for i in range(len(context_information))]
    context_information_metadata_url = [context_information_metadata[i]['url'] for i in range(len(context_information_metadata))]
    context_information_metadata_name = [s[82:] for s in context_information_metadata_url]

    context = str(context_information_metadata_name)
    prompt = output['query']
    if sheet is not None:
        log_to_google_sheet(sheet, context, prompt, prompt_template, response)
    else:
        print("Google Sheet not initialized. Skipping logging...")
    