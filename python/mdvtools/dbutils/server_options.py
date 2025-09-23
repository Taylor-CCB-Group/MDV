from flask import Flask
from mdvtools.server_extension import MDVServerOptions, MDVServerExtension
from mdvtools.llm.chat_server_extension import MDVProjectChatServerExtension
from mdvtools.dbutils.project_manager_extension import ProjectManagerExtension

extensions: dict[str, type[MDVServerExtension]] = {
    "chat": MDVProjectChatServerExtension,
    "project_manager": ProjectManagerExtension,
}


def get_server_options_for_db_projects(app: Flask) -> MDVServerOptions:
    """
    Returns the server options for projects served from the database.
    This includes enabling the chat extension.
    """
    # in future we may have a more structured way of configuring extensions:
    # - bringing in code that isn't part of the mdvtools package
    # - having options to pass to the extensions
    extensions = [extensions[ext]() for ext in app.config['extensions'] if ext in extensions]
    return MDVServerOptions(
        open_browser=False,
        backend_db=True,
        app=app,
        extensions=extensions,
        websocket=True, # to be explicit
    )