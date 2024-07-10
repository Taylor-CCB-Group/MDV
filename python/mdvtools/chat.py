from mdvtools.mdvproject import MDVProject
import os
from dotenv import load_dotenv

load_dotenv()

print(f'Running ChatMDV prototype with OpenAI key {os.getenv("OPENAI_API_KEY")}')

if __name__ == "__main__":
    path = os.path.expanduser("~/mdv/pbmc3k")
    p = MDVProject(path)
    p.serve(websocket=True, use_reloader=True, open_browser=False)
