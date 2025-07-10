from flask import Flask
from mdvtools.server_extension import MDVServerOptions
from mdvtools.llm.chat_server_extension import chat_extension

def get_server_options_for_db_projects(app: Flask) -> MDVServerOptions:
    """
    Returns the server options for projects served from the database.
    This includes enabling the chat extension.
    """
    return MDVServerOptions(
        open_browser=False,
        backend_db=True,
        app=app,
        extensions=[chat_extension],
        websocket=True, # to be explicit
    ) 