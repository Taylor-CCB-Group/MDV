import os
import json
import uuid
import threading
import logging
import tempfile
import time
import shutil
from datetime import datetime, timedelta, timezone
from flask import Flask
from mdvtools.file_processing import datasource_processing
from mdvtools.file_processing import validate_datasource
from flask_sock import Sock
from werkzeug.utils import secure_filename
from typing import Dict, Optional, Callable, Any, List
import traceback
import base64 

# --- Configuration ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

UPLOAD_TTL = timedelta(hours=24)
CLEANUP_INTERVAL_SECONDS = 3600
_mdv_websocket_upload: Optional['MDVWebSocketUpload'] = None

# --- File Upload Manager ---
class FileUploadManager:
    """
    Manages file uploads, including chunked uploads, resumable states,
    finalization, and background cleanup of expired uploads.

    Attributes:
        app (Flask): Flask application context for processing.
        lock (threading.Lock): Lock to control access to shared state.
        processing_queue (List[Dict]): Queue of files waiting to be processed.
        connection_pool (Dict[str, Any]): Mapping of client connections.
        base_temp_dir (str): Directory to store temporary upload files.
        upload_ttl (timedelta): Time to live for upload sessions.
    """
    def __init__(self, app: Flask, base_temp_dir: Optional[str] = None, upload_ttl: timedelta = UPLOAD_TTL):
        """
        Initializes the FileUploadManager with background worker threads.

        Args:
            app (Flask): The Flask application context.
            base_temp_dir (Optional[str], optional): Path for temporary upload storage. Defaults to system temp dir.
            upload_ttl (timedelta, optional): Time to keep uploads before cleanup. Defaults to 24 hours.
        """
        self.app = app
        self.lock = threading.Lock() 
        self.processing_queue: List[Dict] = []
        self.connection_pool: Dict[str, Any] = {} 
        self.base_temp_dir = base_temp_dir or os.path.join(tempfile.gettempdir(), 'mdv_uploads')
        self.upload_ttl = upload_ttl
        os.makedirs(self.base_temp_dir, exist_ok=True)
        logger.info(f"FileUploadManager initialized. Base temp dir: {self.base_temp_dir}, TTL: {self.upload_ttl}")
        self.worker_thread = threading.Thread(target=self._process_queue, daemon=True)
        self.worker_thread.start()
        self.cleanup_stop_event = threading.Event()
        self.cleanup_thread = threading.Thread(target=self._cleanup_abandoned_uploads, daemon=True)
        self.cleanup_thread.start()

    def _get_upload_dir(self, file_id: str) -> str:
        """
        Constructs a secure directory path for storing file uploads.

        Args:
            file_id (str): Unique identifier of the file.

        Returns:
            str: Path to the upload directory.

        Raises:
            ValueError: If the file ID contains invalid characters.
        """
        safe_file_id = secure_filename(file_id)
        if safe_file_id != file_id:
             raise ValueError(f"Invalid file_id format: {file_id}")
        return os.path.join(self.base_temp_dir, safe_file_id)

    def _get_state_filepath(self, file_id: str) -> str:
        """
        Returns the full path to the upload state JSON file.

        Args:
            file_id (str): Unique identifier for the upload session.

        Returns:
            str: Path to the JSON state file.
        """
        return os.path.join(self._get_upload_dir(file_id), 'upload_state.json')

    def _load_upload_state(self, file_id: str) -> Optional[Dict]:
        """
        Loads and parses the upload state for a given file ID.

        Args:
            file_id (str): Unique identifier for the file.

        Returns:
            Optional[Dict]: Upload state dictionary or None on failure.
        """
        state_filepath = self._get_state_filepath(file_id)
        try:
            if os.path.exists(state_filepath):
                with open(state_filepath, 'r') as f:
                    state = json.load(f)
                    if 'last_activity_time' in state and state['last_activity_time'] is not None:
                         state['last_activity_time_iso'] = state['last_activity_time'] 
                         state['last_activity_time'] = datetime.fromisoformat(state['last_activity_time'])
                    return state
        except (json.JSONDecodeError, OSError, TypeError) as e:
            logger.error(f"Error loading state for file_id {file_id}: {e}")
        return None

    def _save_upload_state(self, file_id: str, state: Dict) -> bool:
        """
        Saves the upload state to disk, ensuring data is flushed and synced.

        Args:
            file_id (str): File identifier.
            state (Dict): Upload state data to persist.

        Returns:
            bool: True if saved successfully, False otherwise.
        """
        upload_dir = self._get_upload_dir(file_id)
        state_filepath = self._get_state_filepath(file_id)
        temp_state_filepath = state_filepath + ".tmp"
        try:
            os.makedirs(upload_dir, exist_ok=True) 
            state_copy = state.copy() 
            if 'last_activity_time' in state_copy and isinstance(state_copy['last_activity_time'], datetime):
                 state_copy['last_activity_time'] = state_copy['last_activity_time'].isoformat()
            with open(temp_state_filepath, 'w') as f:
                json.dump(state_copy, f, indent=4)
                f.flush()  
                os.fsync(f.fileno()) 
            os.replace(temp_state_filepath, state_filepath)
            logger.debug(f"Successfully saved and synced state for file_id {file_id} with status {state_copy.get('status')}")
            return True
        except (OSError, TypeError) as e:
            logger.error(f"Error saving state for file_id {file_id}: {e}")
            if os.path.exists(temp_state_filepath):
                try: os.remove(temp_state_filepath)
                except OSError as rm_err: logger.error(f"Error removing temporary state file {temp_state_filepath}: {rm_err}")
            return False

    def _update_last_activity(self, file_id: str, state: Dict) -> None:
        """
        Updates the last activity timestamp in the upload state.

        Args:
            file_id (str): File identifier.
            state (Dict): Current upload state to update.
        """
        state['last_activity_time'] = datetime.now(timezone.utc) 
        self._save_upload_state(file_id, state)

    def _delete_upload_dir(self, file_id: str) -> None:
        """
        Removes the temporary upload directory for a given file ID.

        Args:
            file_id (str): File identifier.
        """
        upload_dir = self._get_upload_dir(file_id)
        if os.path.exists(upload_dir):
            try:
                shutil.rmtree(upload_dir)
                logger.info(f"Cleaned up temporary directory: {upload_dir}")
            except OSError as e:
                logger.error(f"Error deleting directory {upload_dir}: {e}")
        else:
            logger.warning(f"Attempted to delete non-existent directory: {upload_dir}")

    def initialize_upload(self, file_id: str, filename: str, total_size: int, project_id: str, content_type: str, extra_data: Dict) -> Dict:
        """
        Initializes or resumes an upload session, storing the upload state.

        Args:
            file_id (str): Unique identifier for the file.
            filename (str): Original file name.
            total_size (int): Total expected file size in bytes.
            project_id (str): ID of the project receiving the upload.
            content_type (str): MIME type of the file.
            extra_data (Dict): Additional metadata to store with the upload.

        Returns:
            Dict: Upload state dictionary.
        """
        upload_dir = self._get_upload_dir(file_id)
        data_filename = secure_filename(filename) 
        data_filepath = os.path.join(upload_dir, data_filename)
        state = self._load_upload_state(file_id)

        if state: 
            logger.info(f"Resuming/Re-initializing upload for file_id: {file_id}. Original filename: {state.get('original_filename')}, Current: {filename}")
            state['project_id'] = project_id 
            state['filename'] = filename 
            state['total_size'] = total_size 
            state['content_type'] = content_type
            state['data_filepath'] = data_filepath 
            state.update(extra_data) 
            if state.get('status') not in ['completed', 'queued', 'processing']:
                 state['status'] = 'resuming' 
            self._update_last_activity(file_id, state) 
        else: 
            logger.info(f"Initializing new upload for file_id: {file_id}, filename: {filename}")
            os.makedirs(upload_dir, exist_ok=True)
            state = {
                'file_id': file_id, 'filename': filename, 'original_filename': filename, 
                'data_filepath': data_filepath, 'total_size': total_size, 'received_bytes': 0,
                'project_id': project_id, 'content_type': content_type, 'status': 'uploading', 
                'start_time_iso': datetime.now(timezone.utc).isoformat(), **extra_data 
            }
            self._update_last_activity(file_id, state)
        return state

    def write_chunk(self, file_id: str, chunk_data: bytes) -> Optional[Dict]:
        """
        Appends a chunk of data to the upload file and updates the upload state.

        Args:
            file_id (str): ID of the upload session.
            chunk_data (bytes): Chunk of file data to write.

        Returns:
            Optional[Dict]: Updated state if successful, None on error.
        """
        state = self._load_upload_state(file_id)
        if not state or state.get('status') not in ['uploading', 'resuming']:
            logger.warning(f"Received chunk for inactive/invalid upload: {file_id}, status: {state.get('status') if state else 'None'}")
            return None
        try:
            with open(state['data_filepath'], 'ab') as f: f.write(chunk_data)
            state['received_bytes'] += len(chunk_data)
            self._update_last_activity(file_id, state) 
            return state
        except OSError as e:
            logger.error(f"Error writing chunk for file_id {file_id} to {state.get('data_filepath')}: {e}")
            state['status'] = 'error'; state['error_message'] = f"File write error: {e}"
            self._save_upload_state(file_id, state) 
            return None
        except Exception as e: 
             logger.error(f"Unexpected error writing chunk for file_id {file_id}: {e}")
             state['status'] = 'error'; state['error_message'] = f"Unexpected error: {e}"
             self._save_upload_state(file_id, state)
             return None

    def finalize_upload(self, file_id: str) -> Optional[Dict]:
        """
        Finalizes the upload by validating file size and updating status.

        Args:
            file_id (str): ID of the file being finalized.

        Returns:
            Optional[Dict]: Finalized upload state or None on failure.
        """
        state = self._load_upload_state(file_id)
        if not state:
            logger.error(f"Cannot finalize non-existent upload: {file_id}")
            return None
        try:
            actual_size = os.path.getsize(state['data_filepath'])
            expected_size = state['total_size']
            if actual_size != expected_size:
                logger.warning(f"Size mismatch for {file_id}: expected {expected_size}, got {actual_size}. File: {state['data_filepath']}")
                state['status'] = 'error'; state['error_message'] = f"Size mismatch: expected {expected_size}, got {actual_size}"
                self._save_upload_state(file_id, state)
                return None
            else:
                logger.info(f"Upload verified for {file_id}. Size: {actual_size}. File: {state['data_filepath']}")
                state['status'] = 'completed' 
                self._update_last_activity(file_id, state) 
                return state
        except OSError as e:
            logger.error(f"Error verifying file size for {file_id} at {state.get('data_filepath')}: {e}")
            state['status'] = 'error'; state['error_message'] = f"File verification error: {e}"
            self._save_upload_state(file_id, state)
            return None

    def cancel_upload(self, file_id: str) -> None:
        """
        Cancels the upload session and removes its temporary directory.

        Args:
            file_id (str): ID of the upload session to cancel.
        """
        logger.info(f"Cancelling upload and cleaning up for file_id: {file_id}")
        self._delete_upload_dir(file_id)

    def queue_file_for_processing(self, state: Dict, ws: Any) -> None:
        """
        Queues the file for background processing and updates its status.

        Args:
            state (Dict): Upload state of the file.
            ws (Any): WebSocket connection to communicate with the client.
        """
        file_id = state.get('file_id')
        current_status = state.get('status')

        # Allow re-queueing if it was completed, queued, or processing (indicating a potential prior interruption)
        # Or if it's a fresh 'completed' state.
        if current_status not in ['completed', 'queued', 'processing']:
            logger.warning(f"Attempted to queue file {file_id} with invalid status for processing: {current_status}")
            self._notify_client(ws, {'type': 'error', 'file_id': file_id, 'message': f'Cannot queue file with status {current_status}'})
            return

        with self.lock: # Protect access to self.processing_queue
            # Check if already in the current worker's in-memory queue to avoid duplicate processing *by the same worker instance*
            for item in self.processing_queue:
                if item['state']['file_id'] == file_id:
                    logger.info(f"File {file_id} (status: {current_status}) is already in the current processing queue. Will not re-add to this worker's memory queue.")
                    # If it's already in this worker's queue, and its status in state file is 'processing',
                    # it means this worker is already handling it or about to.
                    # If status is 'queued', it means this worker is about to pick it up.
                    # If status is 'completed', it means this worker is about to pick it up and change its status to 'queued'.
                    # No need to change state file here if already in memory queue.
                    return 

            logger.info(f"Queueing file {file_id} (status: {current_status}) for processing.")
            processing_info = { "state": state, "ws": ws }
            self.processing_queue.append(processing_info)
        
        # Update status in the state file to 'queued' to reflect its new state.
        # So if this worker dies, the next worker sees 'queued'.
        if state.get('status') != 'queued': # Avoid redundant save if already 'queued'
            state['status'] = 'queued' 
            self._save_upload_state(file_id, state) 

    def _process_queue(self):
        """
        Background worker method to process uploaded files in a queue.
        Validates project associations and triggers data processing logic.
        """
        global _mdv_websocket_upload 
        while True:
            processing_info = None
            with self.lock:
                if self.processing_queue:
                    processing_info = self.processing_queue.pop(0)
            if processing_info:
                state = processing_info['state']
                ws = processing_info['ws'] 
                file_id = state.get('file_id')
                project_id = state.get('project_id')
                filepath = state.get('data_filepath')
                original_filename = state.get('original_filename') 
                content_type = state.get('content_type')
                
                project = None
                if _mdv_websocket_upload: 
                    project = _mdv_websocket_upload.projects_map.get(project_id)
                
                if not project:
                    logger.error(f"Project instance not found for ID: {project_id} during processing of {file_id}. Cannot process.")
                    state['status'] = 'error'; state['error_message'] = f"Project {project_id} not found for processing."
                    self._save_upload_state(file_id, state)
                    self._notify_client(ws, {'type': 'error', 'file_id': file_id, 'message': state['error_message']})
                    continue 

                logger.info(f"Processing file: {file_id} (type: {content_type}, original: {original_filename}) for project {project_id}")
                state['status'] = 'processing' # Mark as processing before starting
                self._save_upload_state(file_id, state) # Save the 'processing' state
                self._notify_client(ws, {'type': 'processing', 'file_id': file_id, 'message': 'File processing started'})

                try:
                    response_data = {}
                    if content_type == "text/csv":
                        view = state.get('view', 'default')
                        replace = state.get('replace', False)
                        supplied_only = state.get('supplied_only', False) 
                        logger.info(f"Calling datasource_processing for {file_id} with: project_id={project_id}, filepath={filepath}, original_filename={original_filename}, view={view}, replace={replace}, supplied_only={supplied_only}")
                        with self.app.app_context(): 
                            result = datasource_processing(project, filepath, original_filename, view, replace, supplied_only)
                        response_data = {'result': result} 
                    else:
                        logger.warning(f"No specific processor for content type: {content_type} on file {file_id}")
                        response_data = {'warning': f'No processor configured for {content_type}'}
                    
                    # After successful processing, update state to 'processed_success' or similar if needed,
                    # though usually we just delete it.
                    self._notify_client(ws, {
                        'type': 'success', 'file_id': file_id,
                        'message': 'File processed successfully', 'status': 200, **response_data 
                    })
                    logger.info(f"Successfully processed file: {file_id}")
                    self._delete_upload_dir(file_id) # Clean up on success
                except Exception as e:
                    logger.error(f"Error processing file {file_id}: {str(e)}")
                    logger.error(traceback.format_exc())
                    state['status'] = 'error'; state['error_message'] = f"Processing error: {str(e)}"
                    self._save_upload_state(file_id, state) 
                    self._notify_client(ws, {'type': 'error', 'file_id': file_id, 'message': f"Processing error: {str(e)}"})
            else:
                time.sleep(1.0)

    def _cleanup_abandoned_uploads(self):
        """
        Periodically scans and deletes expired or orphaned upload directories
        based on their last activity timestamp or file modification time.
        """
        logger.info("Cleanup thread started.")
        while not self.cleanup_stop_event.is_set():
            try:
                now = datetime.now(timezone.utc)
                for item_name in os.listdir(self.base_temp_dir):
                    item_path = os.path.join(self.base_temp_dir, item_name)
                    if os.path.isdir(item_path):
                        file_id = item_name 
                        state = self._load_upload_state(file_id)
                        if state:
                            status = state.get('status')
                            last_activity_dt = state.get('last_activity_time') 
                            if status in ['uploading', 'resuming', 'error', 'queued', 'processing'] and last_activity_dt: # Added 'processing'
                                if isinstance(last_activity_dt, datetime): 
                                    if now - last_activity_dt > self.upload_ttl:
                                        logger.warning(f"Upload {file_id} (status: {status}) exceeded TTL. Cleaning up.")
                                        self._delete_upload_dir(file_id)
                                else: logger.warning(f"Invalid last_activity_time for {file_id}")
                        else: # Orphaned directory
                            try:
                                dir_mod_time = datetime.fromtimestamp(os.path.getmtime(item_path), timezone.utc)
                                if now - dir_mod_time > self.upload_ttl:
                                     logger.warning(f"Orphaned upload dir: {item_path}. Cleaning up.")
                                     shutil.rmtree(item_path) 
                            except OSError as e: logger.error(f"Error deleting orphaned dir {item_path}: {e}")
            except Exception as e: logger.error(f"Error in cleanup: {e}"); logger.error(traceback.format_exc())
            self.cleanup_stop_event.wait(CLEANUP_INTERVAL_SECONDS)
        logger.info("Cleanup thread stopped.")

    def stop_cleanup(self): 
        """
        Signals the background cleanup thread to stop execution.
        """
        self.cleanup_stop_event.set()

    def _notify_client(self, ws, message: dict) -> None:
        """
        Sends a JSON-encoded message to a WebSocket client.

        Args:
            ws (Any): WebSocket connection object.
            message (dict): Message payload to send.
        """
        if not ws:
            logger.warning(f"WS missing for notify: {message.get('file_id', 'N/A')}")
            return
        try: ws.send(json.dumps(message))
        except Exception as e: logger.error(f"Error sending to client (file_id: {message.get('file_id', 'N/A')}): {e}")

# --- WebSocket Handler ---
class MDVWebSocketUpload:
    """
    WebSocket-based handler for managing real-time file uploads and processing.

    Attributes:
        app (Flask): The Flask application context.
        sock (Sock): WebSocket handler.
        upload_manager (FileUploadManager): Upload manager instance.
        projects_map (Dict[str, Any]): Mapping of project IDs to project instances.
    """
    def __init__(self, app: Flask):
        """
        Initializes the WebSocket handler and binds upload route.

        Args:
            app (Flask): Flask app for WebSocket context.
        """
        self.app = app
        self.sock = Sock(app)
        self.upload_manager = FileUploadManager(app=self.app) 
        self.projects_map: Dict[str, Any] = {} 

        @self.sock.route('/project/<project_id>/ws')
        def handle_websocket_upload(ws, project_id):
            if project_id not in self.projects_map:
                logger.warning(f"WS for unknown project: {project_id}")
                try: ws.send(json.dumps({'type': 'error', 'message': f"Unknown project ID: {project_id}"})); ws.close() 
                except Exception: pass 
                return
            self._handle_project_websocket(ws, project_id)

    def register_project(self, project_id: str, project: Any):
        """
        Registers a new project instance for handling uploads.

        Args:
            project_id (str): Unique identifier of the project.
            project (Any): Project instance to associate.
        """
        self.projects_map[project_id] = project
        logger.info(f"Registered project {project_id} for WebSocket uploads")

    def _handle_project_websocket(self, ws, project_id: str):
        """
        Manages the lifecycle of a WebSocket connection for a specific project,
        handling file upload messages, chunk transmission, finalization, and errors.

        Args:
            ws (Any): WebSocket connection.
            project_id (str): Identifier for the project.
        """
        client_id = str(uuid.uuid4())
        self.upload_manager.connection_pool[client_id] = ws
        current_upload_state: Optional[Dict] = None 
        logger.info(f"WS client connected for project {project_id}: {client_id}")
        self.upload_manager._notify_client(ws, {
            'type': 'connected', 'client_id': client_id,
            'project_id': project_id, 'message': 'Connection established, ready for upload commands.'
        })
        try:
            while True:
                message_str = ws.receive()
                if message_str is None:
                    logger.info(f"WS connection closed by client {client_id} for project {project_id}.")
                    break 
                try:
                    data = json.loads(message_str)
                    msg_type = data.get('type')
                    file_id = data.get('file_id') 

                    if msg_type == 'query_upload':
                        if not file_id: raise ValueError("'query_upload' message requires 'file_id'")
                        logger.info(f"Received query for file_id: {file_id}")
                        state = self.upload_manager._load_upload_state(file_id)
                        if state:
                            status = state.get('status')
                            original_filename = state.get('original_filename', 'unknown')

                            # If file was fully uploaded but potentially interrupted during/before processing
                            if status in ['completed', 'queued', 'processing']:
                                logger.info(f"File {file_id} ({original_filename}) status is '{status}'. Re-initiating processing.")
                                # Re-queue. queue_file_for_processing sets status to 'queued' and saves.
                                # _process_queue will then pick it up, set to 'processing', save, and notify.
                                self.upload_manager.queue_file_for_processing(state, ws) 
                                # Send an immediate notification that reprocessing is being initiated.
                                # The actual 'processing' message will come from _process_queue.
                                self.upload_manager._notify_client(ws, {
                                    'type': 'processing_initiated', # Client should treat this like 'processing'
                                    'file_id': file_id,
                                    'message': f'Processing for {original_filename} is being (re)initiated.'
                                })
                                # Set current_upload_state so that client doesn't try to re-upload chunks
                                current_upload_state = state 
                            elif status in ['uploading', 'resuming', 'error']: 
                                logger.info(f"File {file_id} ({original_filename}) status is '{status}'. Sending resume_info.")
                                self.upload_manager._notify_client(ws, {
                                    'type': 'resume_info', 'file_id': file_id,
                                    'received_bytes': state.get('received_bytes', 0),
                                    'total_size': state.get('total_size', 0) 
                                })
                                current_upload_state = state # Ensure current_upload_state is set for resume
                            else: 
                                 logger.warning(f"File {file_id} ({original_filename}) has unknown status '{status}'. Treating as not found.")
                                 self.upload_manager._notify_client(ws, {'type': 'upload_not_found', 'file_id': file_id})
                        else: 
                            logger.info(f"No state found for file_id {file_id}. Treating as not found.")
                            self.upload_manager._notify_client(ws, {'type': 'upload_not_found', 'file_id': file_id})

                    elif msg_type == 'start':
                        if not file_id:
                            file_id = str(uuid.uuid4())
                            logger.warning(f"Client did not provide file_id for 'start'. Generated new one: {file_id}.")
                        
                        filename = secure_filename(data.get('filename', f'file_{file_id}'))
                        total_size = int(data.get('size', 0))
                        content_type = data.get('content_type', 'application/octet-stream')
                        
                        extra_data = {}
                        project_instance = self.projects_map.get(project_id)
                        if not project_instance: raise ValueError(f"Project instance {project_id} not found.")

                        if content_type == "text/csv":
                             try:
                                 validation_params = {k: data.get(k) for k in ["name", "view", "replace", "supplied_only"]}
                                 validation_params = {k: v for k, v in validation_params.items() if v is not None}
                                 with self.app.app_context():
                                     print(f"Validating datasource with params: {validation_params} - {data}")
                                     validated_data = validate_datasource(project_instance, validation_params)
                                 extra_data.update(validated_data) 
                             except Exception as val_err: 
                                 logger.error(f"Validation failed for CSV upload {file_id}: {val_err}")
                                 self.upload_manager._notify_client(ws, {'type': 'error', 'file_id': file_id, 'message': f"Validation Error: {getattr(val_err, 'message', str(val_err))}"})
                                 continue 
                        
                        current_upload_state = self.upload_manager.initialize_upload(
                            file_id, filename, total_size, project_id, content_type, extra_data
                        )

                        if current_upload_state:
                             ack_type = 'start_ack'
                             message = 'Upload started'
                             # Check if initialize_upload set it to 'resuming' based on loaded state
                             if current_upload_state.get('status') == 'resuming':
                                 ack_type = 'resume_ack' 
                                 message = f"Resuming upload from byte {current_upload_state.get('received_bytes', 0)}"
                             self.upload_manager._notify_client(ws, {
                                'type': ack_type, 'file_id': file_id,
                                'received_bytes': current_upload_state.get('received_bytes', 0), 
                                'message': message
                            })
                             logger.info(f"{message} for file {filename} ({file_id}), size {total_size}, project {project_id}")
                        else:
                             self.upload_manager._notify_client(ws, {'type': 'error', 'file_id': file_id, 'message': 'Failed to initialize upload state.'})

                    elif msg_type == 'chunk':
                        if not current_upload_state or not file_id or current_upload_state.get('file_id') != file_id:
                            logger.warning(f"Chunk for mismatched/missing file_id. Client: {file_id}, Server expects: {current_upload_state.get('file_id') if current_upload_state else 'None'}")
                            self.upload_manager._notify_client(ws, {'type': 'error', 'file_id': file_id, 'message': 'Upload context mismatch. Please restart upload.'})
                            continue

                        chunk_data_b64 = data.get('data')
                        chunk_num = int(data.get('chunk_num', 0)) 
                        try: chunk_data = base64.b64decode(chunk_data_b64)
                        except (TypeError, base64.binascii.Error) as e:
                             logger.error(f"Invalid base64 for chunk {chunk_num} of {file_id}: {e}")
                             self.upload_manager._notify_client(ws, {'type': 'error', 'file_id': file_id, 'message': 'Invalid chunk data encoding.'})
                             if current_upload_state: self.upload_manager.cancel_upload(current_upload_state['file_id'])
                             current_upload_state = None; continue

                        updated_state = self.upload_manager.write_chunk(file_id, chunk_data)
                        if updated_state:
                            current_upload_state = updated_state 
                            if chunk_num % 5 == 0 or current_upload_state['received_bytes'] == current_upload_state['total_size']:
                                progress = min(100, int((current_upload_state['received_bytes'] / current_upload_state['total_size']) * 100)) if current_upload_state['total_size'] > 0 else 100
                                self.upload_manager._notify_client(ws, {
                                    'type': 'progress', 'file_id': file_id,
                                    'received': current_upload_state['received_bytes'],
                                    'total': current_upload_state['total_size'], 'progress': progress
                                })
                        else:
                            logger.error(f"Failed to write chunk {chunk_num} for {file_id}."); current_upload_state = None 

                    elif msg_type == 'end':
                        if not current_upload_state or not file_id or current_upload_state.get('file_id') != file_id:
                            logger.warning(f"END for mismatched/missing file_id. Client: {file_id}, Server: {current_upload_state.get('file_id') if current_upload_state else 'None'}")
                            self.upload_manager._notify_client(ws, {'type': 'error', 'file_id': file_id, 'message': 'Upload context mismatch for end.'})
                            continue
                        final_state = self.upload_manager.finalize_upload(file_id)
                        if final_state:
                            self.upload_manager.queue_file_for_processing(final_state, ws)
                            try:
                                start_time = datetime.fromisoformat(final_state.get('start_time_iso', datetime.now(timezone.utc).isoformat()))
                                elapsed = (datetime.now(timezone.utc) - start_time).total_seconds()
                                speed = final_state['total_size'] / (elapsed or 0.001) / 1024
                            except Exception: elapsed = -1; speed = -1
                            self.upload_manager._notify_client(ws, {
                                'type': 'end_ack', 'file_id': file_id,
                                'message': 'Upload completed and verified. Queued for processing.',
                                'upload_speed_KBps': round(speed, 2) if speed != -1 else 'N/A',
                                'upload_time_sec': round(elapsed, 2) if elapsed != -1 else 'N/A'
                            })
                            logger.info(f"Upload completed and queued: {current_upload_state['filename']} ({file_id})")
                        else: logger.error(f"Upload finalization failed for {file_id}.")
                        current_upload_state = None
                    elif msg_type == 'cancel':
                         if not file_id:
                              logger.warning("Cancel without file_id."); self.upload_manager._notify_client(ws, {'type': 'error', 'message': "'cancel' requires 'file_id'"})
                              continue
                         logger.info(f"Client cancel for file_id: {file_id}"); self.upload_manager.cancel_upload(file_id) 
                         self.upload_manager._notify_client(ws, {'type': 'cancel_ack', 'file_id': file_id, 'message': 'Upload cancelled by client.'})
                         if current_upload_state and current_upload_state.get('file_id') == file_id: current_upload_state = None
                    elif msg_type == 'ping': self.upload_manager._notify_client(ws, {'type': 'pong'})
                    elif msg_type == 'popout': logger.info(f"Popout: {data}"); self.upload_manager._notify_client(ws, {"status": "ok", "action": "popout_received"})
                    else:
                        logger.warning(f"Unknown msg type: {msg_type} for project {project_id}")
                        self.upload_manager._notify_client(ws, {'type': 'error', 'message': f"Unknown message type: {msg_type}"})
                except json.JSONDecodeError:
                    logger.error(f"Invalid JSON from {client_id}: {message_str}")
                    self.upload_manager._notify_client(ws, {'type': 'error', 'message': 'Invalid JSON format'})
                except ValueError as ve: 
                     logger.error(f"Value error from {client_id} (file_id: {file_id if file_id else 'N/A'}): {ve}")
                     self.upload_manager._notify_client(ws, {'type': 'error', 'file_id': file_id, 'message': f"Invalid request: {ve}"})
                except Exception as e:
                    logger.exception(f"Unexpected error handling WS msg from {client_id} for project {project_id} (file_id: {file_id if file_id else 'N/A'})")
                    self.upload_manager._notify_client(ws, {'type': 'error', 'file_id': file_id, 'message': 'Internal server error.'})
        except Exception as e:
            if "connection closed" not in str(e).lower() and "socket closed" not in str(e).lower():
                logger.exception(f"WS connection error for client {client_id}, project {project_id}: {e}")
            else: logger.info(f"WS connection gracefully closed for client {client_id}, project {project_id}.")
        finally:
            if client_id in self.upload_manager.connection_pool: del self.upload_manager.connection_pool[client_id]
            logger.info(f"Finished handling WS connection for client {client_id}, project {project_id}")

# --- Singleton and Initialization ---
_mdv_websocket_upload: Optional[MDVWebSocketUpload] = None

def initialize_websocket_upload(app: Flask) -> MDVWebSocketUpload:
    """
    Initializes the singleton instance of MDVWebSocketUpload if not already created.

    Args:
        app (Flask): The Flask application instance.

    Returns:
        MDVWebSocketUpload: The initialized WebSocket upload handler instance.
    """
    global _mdv_websocket_upload
    if _mdv_websocket_upload is None:
        _mdv_websocket_upload = MDVWebSocketUpload(app)
        logger.info("MDVWebSocketUpload initialized successfully.")
    else: logger.warning("MDVWebSocketUpload already initialized.")
    return _mdv_websocket_upload

def register_project_for_upload(project_id: str, project: Any):
    """
    Registers a project to be available for file uploads via WebSocket.

    Args:
        project_id (str): Unique identifier of the project.
        project (Any): Project instance to associate with the given ID.

    Raises:
        RuntimeError: If the WebSocket upload manager has not been initialized.
    """
    global _mdv_websocket_upload
    if _mdv_websocket_upload is None:
        raise RuntimeError("WebSocket upload manager not initialized. Call initialize_websocket_upload first.")
    _mdv_websocket_upload.register_project(project_id, project)
