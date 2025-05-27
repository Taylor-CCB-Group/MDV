from flask import Flask, request
from flask_socketio import SocketIO
from datetime import datetime
# import asyncio
# from typing import Optional
from mdvtools.mdvproject import MDVProject

import logging

def log(msg: str):
    now = datetime.now()
    date_str = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[[[ socket.io ]]] [{date_str}] - {msg}")

socketio: SocketIO = None # type: ignore

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
    cors = [f"http://localhost:{i}" for i in range(5170,5180)]
    socketio = SocketIO(app, cors_allowed_origins=cors)
    log("socketio initialized without cors_allowed_origins wildcard")

    async def response(sid, message=""):
        # await asyncio.sleep(1)
        socketio.emit("message", {"type": "ping", "message": f"bleep bloop {message} I'm a robot"}, to=sid)

    @socketio.on("connect")
    def test_connect(auth):
        log(f"WebSocket client connected, args: {request.args}")
        # socketio.emit('popout', "KNkyOT")

    @socketio.on("message")
    def message(data):
        # 'hello world'... if we have need for more logic here, we'll probably use a dictionary of functions for each message type.
        if data['type'] == "popout":
            log(f"popout: {data}")
        if data['type'] == "ping":
            log(f"ping: {data}")

            # socketio.emit("message", {"type": "ping", "message": "bleep bloop I'm a robot"})
            # how about getting it to respond only to the client that sent the message?
            # socketio.emit("message", {"type": "ping", "message": "bleep bloop I'm a robot"}, to=request.sid)
            # ... or maybe I should stick to REST for now.
            sid = request.args.get("sid")
            # error asyncio.run() cannot be called from a running event loop
            # asyncio.run(response(sid, data.get('message')))
            # response(sid, data.get('message'))


def get_socket_logger(event_name="log"):
    logger = logging.getLogger(event_name)
    logger.setLevel(logging.INFO)
    handler = SocketIOHandler(socketio, event_name)
    logger.addHandler(handler)
    return logger

class ChatSocketAPI:
    def __init__(self, project: MDVProject):
        self.project = project
        self.socketio = socketio
        # todo refactor event/room/namespace names 
        log_name = f"/project/{project.id}/chat"
        self.progress_name = f"/project/{project.id}/chat_progress"
        self.logger = get_socket_logger(log_name)
    def update_chat_progress(self, message: str, id: str, progress: int, delta: int):
        """
        Send a message to the chat log and also update the progress bar.
        
        Args:
            message (str): the message to send to the chat log
            id (str): the id of the associated chat request
            progress (int): the progress value (%) to update the progress bar with
            delta (int): the expected cost of the current operation (%)
        """
        # we should descriminate which user to send this to... 
        # which implies that this instance should be associated...
        self.socketio.emit(self.progress_name, {"message": message, "id": id, "progress": progress, "delta": delta})

class SocketIOHandler(logging.StreamHandler):
    def __init__(self, socketio, event_name="log"):
        super().__init__()
        log(f"handler initialized for event: {event_name}")
        self.socketio = socketio
        # todo - event_name vs namespace vs room refactor
        self.event_name = event_name
        def my_function_handler(data):
            log(f"{event_name}: {data}")
        # log(f"socketio.on_event({event_name}, my_function_handler)")
        socketio.on_event(event_name, my_function_handler)
        # socketio.on(event_name, )

    def emit(self, record):
        """
        Emit a record - send it via socketio & also print it to the console.
        subject to change.
        """
        try:
            msg = self.format(record)
            log(f"[ {self.event_name} ] {msg}")
            self.socketio.emit(self.event_name, msg)
        except Exception:
            self.handleError(record)
