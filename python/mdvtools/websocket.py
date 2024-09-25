from flask import Flask
from flask_socketio import SocketIO
from datetime import datetime


def log(msg: str):
    now = datetime.now()
    date_str = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[[[ socket.io ]]] [{date_str}] - {msg}")


def mdv_socketio(app: Flask):
    socketio = SocketIO(app, cors_allowed_origins="*")

    @socketio.on("connect")
    def test_connect():
        log("WebSocket client connected")
        # socketio.emit('popout', "KNkyOT")

    @socketio.on("message")
    def message(data):
        # 'hello world'... if we have need for more logic here, we'll probably use a dictionay of functions for each message type.
        if data.type == "popout":
            log(f"popout: {data}")
