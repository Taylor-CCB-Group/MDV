from flask import Flask
from mdvtools.server_extension import MDVServerOptions, MDVProjectServerExtension
from mdvtools.llm.chat_server_extension import MDVProjectChatServerExtension
from mdvtools.dbutils.project_manager_extension import ProjectManagerExtension
from mdvtools.ucsc_proxy_extension import UcscProxyServerExtension

extension_classes: dict[str, type[MDVProjectServerExtension]] = {
    "chat": MDVProjectChatServerExtension,
    "project_manager": ProjectManagerExtension,
    # app-wide integrations that are safe to enable by default
    "ucsc_proxy": UcscProxyServerExtension,
}


def get_server_options_for_db_projects(app: Flask) -> MDVServerOptions:
    """
    Returns the server options for projects served from the database.
    This  includes enabling configured extensions like chat and project_manager.
    """
    # in future we may have a more structured way of configuring extensions:
    # - bringing in code that isn't part of the mdvtools package
    # - having options to pass to the extensions
    # Always register the UCSC proxy so the frontend works consistently across:
    # - db-backed deployments (multi-project)
    # - single-project servers
    extensions = [UcscProxyServerExtension()]
    for ext_name in app.config.get('extensions', []):
        if ext_name in extension_classes:
            # Avoid double-instantiating if the extension is explicitly enabled.
            if ext_name != "ucsc_proxy":
                extensions.append(extension_classes[ext_name]())
        else:
            # Log warning for unknown extensions but don't fail
            import logging
            logging.warning(f"Unknown extension '{ext_name}' in config, skipping")
    
    return MDVServerOptions(
        open_browser=False,
        backend_db=True,
        app=app,
        extensions=extensions,
        websocket=True, # to be explicit
    )