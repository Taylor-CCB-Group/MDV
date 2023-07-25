from flask import Flask
from flask_socketio import SocketIO
from datetime import datetime

def log(msg: str):
    now = datetime.now()
    date_str = now.strftime('%Y-%m-%d %H:%M:%S')
    print(f'[[[ socket.io ]]] [{date_str}] - {msg}')

def mdv_socketio(app: Flask):
    socketio = SocketIO(app, cors_allowed_origins="*")
    @socketio.on('connect')
    def test_connect():
        log('WebSocket client connected')
        # socketio.emit('popout', "KNkyOT")
    
    @socketio.on('popout')
    def popout(data):
        log(f'popout: {data}')