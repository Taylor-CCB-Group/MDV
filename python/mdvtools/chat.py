from mdvtools.mdvproject import MDVProject
import os
from dotenv import load_dotenv
from langchain_openai import ChatOpenAI
import json
# make sure there's a .env file in the same directory as this script with OPENAI_API_KEY=... & GITHUB_TOKEN=...
load_dotenv()

code_llm = ChatOpenAI(temperature=0.1, model_name="gpt-4o")
dataframe_llm = ChatOpenAI(temperature=0.1, model_name="gpt-4")

def chat_fn(question: str) -> str:
    response = code_llm.invoke(question)
    return json.dumps(response.to_json())


if __name__ == "__main__":
    path = os.path.expanduser("~/mdv/pbmc3k")
    p = MDVProject(path)
    p.serve(websocket=True, use_reloader=True, open_browser=False, chat_fn=chat_fn)
