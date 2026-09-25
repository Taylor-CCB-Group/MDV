from flask import Flask, Request as FlaskRequest
from flask_socketio import SocketIO
from datetime import datetime
from typing import Optional

def log(msg: str):
    now = datetime.now()
    date_str = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[[[ socket.io ]]] [{date_str}] - {msg}")

socketio: Optional[SocketIO] = None

# This is a map of chat request IDs to user IDs, used to track which user is associated with which chat session.
# something of a placeholder
# chat_sid_map: dict[str, str] = {}


class SocketIOContextRequest(FlaskRequest):
    """
    Represents the Flask request object when in a Flask-SocketIO event context.
    This protocol ensures the 'sid' (session ID) attribute is recognized by type checkers.
    """
    sid: str
    # namespace: Optional[str]
    # event: Optional[str]

def mdv_socketio(app: Flask):
    """
    Experimental and not to be trusted pending design work etc.

    Do we have a SocketIO for the entire app and route messages internally,
    - yes probably.
    What `path` should it use? We'll need to make sure it is compatible with
    - the vite dev server proxy
    - the app not being at the root of the domain

    !! What about cors? We had a wildcard for previous experiment - should be reviewed before getting pushed online.
    ^^^^
    """
    global socketio
    # allow cors for localhost:5170-5179
    # cors = [f"http://localhost:{i}" for i in range(5170,5180)]

    # Get the API root path for SocketIO path configuration
    api_root = app.config.get('mdv_api_root', '/')
    if api_root != '/' and not api_root.endswith('/'):
        api_root += '/'
    socketio_path = f"{api_root}socket.io"

    socketio = SocketIO(app, cors_allowed_origins="*", path=socketio_path)
    log(f"socketio initialized with cors_allowed_origins wildcard and path={socketio_path}")

    # @socketio.on("message")
    # def message(data):
    #     # 'hello world'... if we have need for more logic here, we'll probably use a dictionary of functions for each message type.
    
    #     if data['type'] == "popout":
    #         log(f"popout: {data}")
    #     if data['type'] == "ping":
    #         log(f"ping: {data}")

    #         # socketio.emit("message", {"type": "ping", "message": "bleep bloop I'm a robot"})
    #         # how about getting it to respond only to the client that sent the message?
    #         # socketio.emit("message", {"type": "ping", "message": "bleep bloop I'm a robot"}, to=request.sid)
    #         # ... or maybe I should stick to REST for now.
    #         sid = request.args.get("sid")
    #         # error asyncio.run() cannot be called from a running event loop
    #         # asyncio.run(response(sid, data.get('message')))
    #         # response(sid, data.get('message'))
