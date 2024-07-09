from flask import Flask, request
from flask_socketio import SocketIO
from datetime import datetime
import asyncio

def log(msg: str):
    now = datetime.now()
    date_str = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[[[ socket.io ]]] [{date_str}] - {msg}")


def mdv_socketio(app: Flask):
    socketio = SocketIO(app, cors_allowed_origins="*")
    print("socketio initialized")

    async def response(sid, message=""):
        await asyncio.sleep(1)
        socketio.emit("message", {"type": "ping", "message": f"bleep bloop {message} I'm a robot"}, to=sid)

    @socketio.on("connect")
    def test_connect(auth):
        log(f"WebSocket client connected, args: {request.args}")
        # socketio.emit('popout', "KNkyOT")

    @socketio.on("message")
    def message(data):
        # 'hello world'... if we have need for more logic here, we'll probably use a dictionay of functions for each message type.
        if data['type'] == "popout":
            log(f"popout: {data}")
        if data['type'] == "ping":
            log(f"ping: {data}")

            # socketio.emit("message", {"type": "ping", "message": "bleep bloop I'm a robot"})
            # how about getting it to respond only to the client that sent the message?
            # socketio.emit("message", {"type": "ping", "message": "bleep bloop I'm a robot"}, to=request.sid)
            # ... or maybe I should stick to REST for now.
            sid = request.args.get("sid")
            asyncio.run(response(sid, data.get('message')))

