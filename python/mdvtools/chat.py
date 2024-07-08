from mdvtools.mdvproject import MDVProject
import os


if __name__ == "__main__":
    path = os.path.expanduser("~/mdv/pbmc3k")
    p = MDVProject(path)
    p.serve(websocket=True, use_reloader=True, open_browser=False)
