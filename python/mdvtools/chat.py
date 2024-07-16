from mdvtools.mdvproject import MDVProject
import os
from dotenv import load_dotenv
from langchain_openai import ChatOpenAI
# make sure there's a .env file in the same directory as this script with OPENAI_API_KEY=... & GITHUB_TOKEN=...
load_dotenv()

# consider what context the agents will run in... for now, this global thing that is only running on a given test project...
# but we'd expect it to be in the context of a particular session interacting with a particular project
code_llm = ChatOpenAI(temperature=0.1, model_name="gpt-4o")
dataframe_llm = ChatOpenAI(temperature=0.1, model_name="gpt-4")

def chat_fn(question: str) -> str:
    response = code_llm.invoke(question)
    # return json.dumps(response.to_json())
    return response.content

welcome = "Hello! I'm a generic gpt-4o chatbot. Ask me anything!"

if __name__ == "__main__":
    path = os.path.expanduser("~/mdv/pbmc3k")
    p = MDVProject(path)
    p.serve(websocket=True, use_reloader=True, open_browser=False, chat_fn=chat_fn, chat_welcome=welcome)
