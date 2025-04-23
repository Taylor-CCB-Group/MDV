from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time
import threading
import os

class ProjectDirectoryHandler(FileSystemEventHandler):
    def __init__(self, app, watch_path):
        self.app = app
        self.watch_path = watch_path

    def on_created(self, event):
        if event.is_directory:
            print(f"[Watcher] New project directory detected: {event.src_path}")
            with self.app.app_context():
                serve_projects_from_filesystem(self.app, self.watch_path)

    def on_moved(self, event):
        if event.is_directory:
            print(f"[Watcher] Project directory moved: {event.dest_path}")
            with self.app.app_context():
                serve_projects_from_filesystem(self.app, self.watch_path)

def start_watching(app, path):
    print(f"[Watcher] Starting to watch {path}")
    event_handler = ProjectDirectoryHandler(app, path)
    observer = Observer()
    observer.schedule(event_handler, path=path, recursive=False)
    observer.daemon = True
    observer.start()
