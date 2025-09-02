import os
import json
import uuid
import threading
import tempfile
import time
import shutil
import base64
from datetime import datetime, timedelta, timezone
from flask import Flask, request
from flask_socketio import SocketIO, emit
from werkzeug.utils import secure_filename
from typing import Dict, Optional, Any, List

from mdvtools.mdvproject import MDVProject
from mdvtools.file_processing import datasource_processing, validate_datasource, anndata_processing, validate_anndata, mdv_project_processing, validate_mdv_project

# Upload configuration
UPLOAD_TTL = timedelta(hours=24)
CLEANUP_INTERVAL_SECONDS = 3600

# Global state
_upload_manager: Optional['FileUploadManager'] = None
_upload_projects_map: Dict[str, MDVProject] = {}
_socketio_instance: Optional[SocketIO] = None
_project_creation_api: Optional['ProjectCreationUploadAPI'] = None

def upload_log(msg: str):
    """Logging function for upload module."""
    now = datetime.now()
    date_str = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[[[ socketio-upload ]]] [{date_str}] - {msg}")

class FileUploadManager:
    """
    Manages file uploads, including chunked uploads, resumable states,
    finalization, and background cleanup of expired uploads.
    """
    
    def __init__(self, app: Flask, socketio: SocketIO, base_temp_dir: Optional[str] = None, upload_ttl: timedelta = UPLOAD_TTL):
        """Initialize the FileUploadManager with background worker threads."""
        self.app = app
        self.socketio = socketio
        self.lock = threading.Lock() 
        self.processing_queue: List[Dict] = []
        self.base_temp_dir = base_temp_dir or os.path.join(tempfile.gettempdir(), 'mdv_uploads')
        self.upload_ttl = upload_ttl
        os.makedirs(self.base_temp_dir, exist_ok=True)
        upload_log(f"FileUploadManager initialized. Base temp dir: {self.base_temp_dir}, TTL: {self.upload_ttl}")
        
        # Start background workers
        self.worker_thread = threading.Thread(target=self._process_queue, daemon=True)
        self.worker_thread.start()
        self.cleanup_stop_event = threading.Event()
        self.cleanup_thread = threading.Thread(target=self._cleanup_abandoned_uploads, daemon=True)
        self.cleanup_thread.start()

    def _get_upload_dir(self, file_id: str) -> str:
        """Constructs a secure directory path for storing file uploads."""
        safe_file_id = secure_filename(file_id)
        if safe_file_id != file_id:
            raise ValueError(f"Invalid file_id format: {file_id}")
        return os.path.join(self.base_temp_dir, safe_file_id)

    def _get_state_filepath(self, file_id: str) -> str:
        """Returns the full path to the upload state JSON file."""
        return os.path.join(self._get_upload_dir(file_id), 'upload_state.json')

    def _load_upload_state(self, file_id: str) -> Optional[Dict]:
        """Loads and parses the upload state for a given file ID."""
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
            upload_log(f"Error loading state for file_id {file_id}: {e}")
        return None

    def _save_upload_state(self, file_id: str, state: Dict) -> bool:
        """Saves the upload state to disk, ensuring data is flushed and synced."""
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
            upload_log(f"Successfully saved state for file_id {file_id} with status {state_copy.get('status')}")
            return True
        except (OSError, TypeError) as e:
            upload_log(f"Error saving state for file_id {file_id}: {e}")
            if os.path.exists(temp_state_filepath):
                try: os.remove(temp_state_filepath)
                except OSError: pass
            return False

    def _update_last_activity(self, file_id: str, state: Dict) -> None:
        """Updates the last activity timestamp in the upload state."""
        state['last_activity_time'] = datetime.now(timezone.utc) 
        self._save_upload_state(file_id, state)

    def _delete_upload_dir(self, file_id: str) -> None:
        """Removes the temporary upload directory for a given file ID."""
        upload_dir = self._get_upload_dir(file_id)
        if os.path.exists(upload_dir):
            try:
                shutil.rmtree(upload_dir)
                upload_log(f"Cleaned up temporary directory: {upload_dir}")
            except OSError as e:
                upload_log(f"Error deleting directory {upload_dir}: {e}")

    def initialize_upload(self, file_id: str, filename: str, total_size: int, project_id: str, content_type: str, extra_data: Dict) -> Dict:
        """Initialize or resume an upload session."""
        upload_dir = self._get_upload_dir(file_id)
        data_filename = secure_filename(filename)
        data_filepath = os.path.join(upload_dir, data_filename)
        state = self._load_upload_state(file_id)

        if state:
            upload_log(f"Resuming/Re-initializing upload for file_id: {file_id}")
            # Update metadata from the new 'start' message, but preserve key resume state
            state['project_id'] = project_id
            state['filename'] = filename
            state['total_size'] = total_size
            state['content_type'] = content_type
            state['data_filepath'] = data_filepath
            state.update(extra_data)

            # Only adjust file and received_bytes if not already completed/queued/processing
            if state.get('status') not in ['completed', 'queued', 'processing']:
                state['status'] = 'resuming'
                expected_offset = state.get('received_bytes', 0)

                if os.path.exists(data_filepath):
                    try:
                        actual_size = os.path.getsize(data_filepath)
                        if actual_size > expected_offset:
                            upload_log(f"Resuming {file_id}: Truncating file from {actual_size} to {expected_offset}")
                            with open(data_filepath, 'r+b') as f:
                                f.truncate(expected_offset)
                        elif actual_size < expected_offset:
                            upload_log(f"Resuming {file_id}: File size inconsistent. Resetting.")
                            state['received_bytes'] = 0
                            try:
                                os.remove(data_filepath)
                            except OSError:
                                pass
                    except OSError as e:
                        upload_log(f"Error during resume file check: {e}")
                        state['status'] = 'error'
                        state['error_message'] = f"File system error on resume: {e}"

            self._update_last_activity(file_id, state)
        else:
            upload_log(f"Initializing new upload for file_id: {file_id}, filename: {filename}")
            os.makedirs(upload_dir, exist_ok=True)

            # Clean up any orphaned files
            if os.path.exists(data_filepath):
                try:
                    os.remove(data_filepath)
                except OSError:
                    pass

            state = {
                'file_id': file_id, 'filename': filename, 'original_filename': filename,
                'data_filepath': data_filepath, 'total_size': total_size, 'received_bytes': 0,
                'project_id': project_id, 'content_type': content_type, 'status': 'uploading',
                'start_time_iso': datetime.now(timezone.utc).isoformat(), **extra_data
            }
            self._update_last_activity(file_id, state)
        return state

    def write_chunk(self, file_id: str, chunk_data: bytes) -> Optional[Dict]:
        """Write a chunk of data to the upload file."""
        state = self._load_upload_state(file_id)
        if not state or state.get('status') not in ['uploading', 'resuming']:
            upload_log(f"Invalid state for chunk write: {file_id}")
            return None

        if state['received_bytes'] >= state['total_size']:
            upload_log(f"File {file_id}: Already complete, ignoring chunk")
            self._update_last_activity(file_id, state)
            return state

        try:
            with open(state['data_filepath'], 'ab') as f:
                f.write(chunk_data)
            
            state['received_bytes'] += len(chunk_data)
            self._update_last_activity(file_id, state)
            return state
        except OSError as e:
            upload_log(f"Error writing chunk for {file_id}: {e}")
            state['status'] = 'error'
            state['error_message'] = f"File write error: {e}"
            self._save_upload_state(file_id, state)
            return None

    def finalize_upload(self, file_id: str) -> Optional[Dict]:
        """Finalize the upload by validating file size."""
        state = self._load_upload_state(file_id)
        if not state:
            return None
            
        try:
            actual_size = os.path.getsize(state['data_filepath'])
            expected_size = state['total_size']

            if actual_size > expected_size:
                upload_log(f"Truncating oversized file {file_id}: {actual_size} -> {expected_size}")
                with open(state['data_filepath'], 'r+b') as f:
                    f.truncate(expected_size)
                actual_size = os.path.getsize(state['data_filepath'])

            if actual_size != expected_size:
                upload_log(f"Size mismatch for {file_id}: expected {expected_size}, got {actual_size}")
                state['status'] = 'error'
                state['error_message'] = f"Size mismatch: expected {expected_size}, got {actual_size}"
                self._save_upload_state(file_id, state)
                return None
            else:
                state['status'] = 'completed'
                self._update_last_activity(file_id, state)
                return state
        except OSError as e:
            upload_log(f"Error verifying file size for {file_id}: {e}")
            state['status'] = 'error'
            state['error_message'] = f"File verification error: {e}"
            self._save_upload_state(file_id, state)
            return None

    def queue_file_for_processing(self, state: Dict, sid: Optional[str], namespace: str) -> None:
        """Queue the file for background processing."""
        file_id = state.get('file_id')
        if not file_id:
            upload_log("Cannot queue file: file_id is missing from state")
            return
        current_status = state.get('status')

        if current_status not in ['completed', 'queued', 'processing']:
            upload_log(f"Cannot queue file {file_id} with status {current_status}")
            self.socketio.emit('upload_error', {
                'file_id': file_id, 'message': f'Cannot queue file with status {current_status}'
            }, namespace=namespace)
            return

        with self.lock:
            # Check if already in queue
            for item in self.processing_queue:
                if item['state']['file_id'] == file_id:
                    upload_log(f"File {file_id} already in processing queue")
                    return 

            upload_log(f"Queueing file {file_id} for processing")
            processing_info = {"state": state, "sid": sid, "namespace": namespace}
            self.processing_queue.append(processing_info)
        
        if state.get('status') != 'queued':
            state['status'] = 'queued'
            if file_id:
                self._save_upload_state(file_id, state)

    def _process_queue(self):
        """Background worker to process uploaded files."""
        while True:
            processing_info = None
            with self.lock:
                if self.processing_queue:
                    processing_info = self.processing_queue.pop(0)
                    
            if processing_info:
                state = processing_info['state']
                sid = processing_info['sid']
                namespace = processing_info['namespace']
                file_id = state.get('file_id')
                # This should never happen if queue_file_for_processing validates properly
                assert file_id, f"file_id missing from queued state: {state}"
                project_id = state.get('project_id')
                
                # Get project from global projects map
                project = _upload_projects_map.get(project_id)
                
                if not project:
                    upload_log(f"Project {project_id} not found for processing {file_id}")
                    
                    # Check if this is a project creation upload
                    if project_id == "project_creation":
                        upload_log(f"Processing project creation upload: {file_id}")
                        # This is handled in the processing logic below
                    else:
                        # file_id is guaranteed to exist due to assertion above
                        state['status'] = 'error'
                        state['error_message'] = f"Project {project_id} not found"
                        self._save_upload_state(file_id, state)
                        self.socketio.emit('upload_error', {
                            'file_id': file_id, 'message': state['error_message']
                        }, namespace=namespace)
                        continue

                upload_log(f"Processing file: {file_id} for project {project_id}")
                state['status'] = 'processing'
                # file_id is guaranteed to exist due to assertion above
                self._save_upload_state(file_id, state)
                self.socketio.emit('upload_processing', {
                    'file_id': file_id, 'message': 'File processing started'
                }, namespace=namespace)

                try:
                    response_data = {}
                    content_type = state.get('content_type')
                    
                    if content_type == "text/csv":
                        # Get project from global projects map
                        project = _upload_projects_map.get(project_id)
                        if not project:
                            raise Exception(f"Project {project_id} not found")
                        
                        # Copy the CSV file to the project directory before processing
                        project_csv_path = os.path.join(project.dir, state['original_filename'])
                        shutil.copy2(state['data_filepath'], project_csv_path)
                        upload_log(f"Copied CSV file to project directory: {project_csv_path}")
                            
                        with self.app.app_context():
                            result = datasource_processing(
                                project, 
                                state['data_filepath'], 
                                state['original_filename'],
                                state.get('view', 'default'),
                                state.get('replace', False),
                                state.get('supplied_only', False)
                            )
                        response_data = {'result': result}
                    elif content_type == "application/x-hdf" or state['original_filename'].endswith('.h5ad'):
                        # Process AnnData file
                        project = _upload_projects_map.get(project_id)
                        if not project:
                            raise Exception(f"Project {project_id} not found")
                            
                        # Copy the h5ad file to the project directory before processing
                        project_h5ad_path = os.path.join(project.dir, state['original_filename'])
                        shutil.copy2(state['data_filepath'], project_h5ad_path)
                        upload_log(f"Copied h5ad file to project directory: {project_h5ad_path}")
                            
                        with self.app.app_context():
                            result = anndata_processing(
                                project,
                                state['data_filepath'],
                                state['original_filename']
                            )
                        response_data = {'result': result}
                    elif (content_type == "application/zip" or state['original_filename'].endswith('.zip')) and project_id == "project_creation":
                        # Process MDV project zip file for project creation
                        with self.app.app_context():
                            # Get projects base directory from app config
                            projects_base_dir = self.app.config.get('projects_base_dir')
                            if not projects_base_dir:
                                raise Exception("Projects base directory not configured")
                            
                            result = mdv_project_processing(
                                self.app,  # Pass the Flask app instance
                                projects_base_dir,
                                state['data_filepath'],
                                state['original_filename'],
                                state.get('project_name')
                            )
                        response_data = {'result': result}
                    else:
                        response_data = {'warning': f'No processor for {content_type} with project_id {project_id}'}
                    
                    self.socketio.emit('upload_success', {
                        'file_id': file_id,
                        'message': 'File processed successfully',
                        **response_data 
                    }, namespace=namespace)
                    upload_log(f"Successfully processed file: {file_id}")
                    self._delete_upload_dir(file_id)
                    
                except Exception as e:
                    upload_log(f"Error processing file {file_id}: {str(e)}")
                    state['status'] = 'error'
                    state['error_message'] = f"Processing error: {str(e)}"
                    # file_id is guaranteed to exist due to assertion above
                    self._save_upload_state(file_id, state)
                    try:
                        self.socketio.emit('upload_error', {
                            'file_id': file_id, 'message': f"Processing error: {str(e)}"
                        }, namespace=namespace)
                    except Exception as emit_error:
                        upload_log(f"Failed to emit error for {file_id}: {emit_error}")
            else:
                time.sleep(1.0)

    def _cleanup_abandoned_uploads(self):
        """Periodically clean up expired uploads."""
        upload_log("Upload cleanup thread started")
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
                            if status in ['uploading', 'resuming', 'error', 'queued', 'processing'] and last_activity_dt:
                                if isinstance(last_activity_dt, datetime): 
                                    if now - last_activity_dt > self.upload_ttl:
                                        upload_log(f"Upload {file_id} expired, cleaning up")
                                        self._delete_upload_dir(file_id)
                        else:
                            # Orphaned directory
                            try:
                                dir_mod_time = datetime.fromtimestamp(os.path.getmtime(item_path), timezone.utc)
                                if now - dir_mod_time > self.upload_ttl:
                                    upload_log(f"Cleaning up orphaned directory: {item_path}")
                                    shutil.rmtree(item_path) 
                            except OSError:
                                pass
            except Exception as e:
                upload_log(f"Error in cleanup: {e}")
            self.cleanup_stop_event.wait(CLEANUP_INTERVAL_SECONDS)

    def cancel_upload(self, file_id: str) -> None:
        """Cancel upload and clean up."""
        upload_log(f"Cancelling upload: {file_id}")
        self._delete_upload_dir(file_id)

    def _notify_client(self, sid: str, namespace: str, event: str, data: dict) -> None:
        """Send a message to a specific client via SocketIO."""
        try:
            self.socketio.emit(event, data, to=sid, namespace=namespace)
        except Exception as e:
            upload_log(f"Error notifying client {sid}: {e}")

class ProjectCreationUploadAPI:
    """SocketIO-based file upload API for creating new projects (zip uploads)."""
    
    def __init__(self, app: Flask, socketio: SocketIO, upload_manager: FileUploadManager):
        self.app = app
        self.socketio = socketio
        self.upload_manager = upload_manager
        self.namespace = "/"  # Use main namespace for project creation
        
        # Register event handlers
        self._register_upload_events()
        upload_log(f"ProjectCreationUploadAPI initialized for project creation uploads")

    def _register_upload_events(self):
        """Register SocketIO event handlers for project creation upload operations."""
        
        @self.socketio.on('upload_query', namespace=self.namespace)
        def handle_upload_query(data):
            """Handle upload status query for project creation."""
            file_id = data.get('file_id')
            if not file_id:
                emit('upload_error', {'message': 'file_id required for query'})
                return
                
            upload_log(f"Project creation upload query for file_id: {file_id}")
            state = self.upload_manager._load_upload_state(file_id)
            
            if state:
                status = state.get('status')
                # Only handle zip files for project creation
                if not (state.get('content_type') == "application/zip" or 
                        state.get('original_filename', '').endswith('.zip')):
                    emit('upload_error', {'message': 'This namespace only handles project zip uploads'})
                    return
                    
                if status in ['completed', 'queued', 'processing']:
                    upload_log(f"File {file_id} status: {status}, re-initiating processing")
                    self.upload_manager.queue_file_for_processing(state, None, self.namespace)
                    emit('upload_processing_initiated', {
                        'file_id': file_id,
                        'message': f'Processing for {state.get("original_filename", "file")} initiated'
                    })
                elif status in ['uploading', 'resuming', 'error']:
                    emit('upload_resume_info', {
                        'file_id': file_id,
                        'received_bytes': state.get('received_bytes', 0),
                        'total_size': state.get('total_size', 0)
                    })
                else:
                    emit('upload_not_found', {'file_id': file_id})
            else:
                emit('upload_not_found', {'file_id': file_id})

        @self.socketio.on('upload_start', namespace=self.namespace)
        def handle_upload_start(data):
            """Handle upload start for project creation."""
            try:
                file_id = data.get('file_id') or str(uuid.uuid4())
                filename = secure_filename(data.get('filename', f'file_{file_id}'))
                total_size = int(data.get('size', 0))
                content_type = data.get('content_type', 'application/octet-stream')
                
                # Only allow zip files for project creation
                if not (content_type == "application/zip" or filename.endswith('.zip')):
                    emit('upload_error', {
                        'file_id': file_id,
                        'message': 'This namespace only accepts zip files for project creation'
                    })
                    return
                
                extra_data = {}
                
                # Validate MDV project zip file
                try:
                    validation_params = {k: data.get(k) for k in ["project_name"]}
                    validation_params = {k: v for k, v in validation_params.items() if v is not None}
                    with self.upload_manager.app.app_context():
                        projects_base_dir = self.upload_manager.app.config.get('projects_base_dir')
                        if not projects_base_dir:
                            raise Exception("Projects base directory not configured")
                        validated_data = validate_mdv_project(self.upload_manager.app, projects_base_dir, validation_params)
                    extra_data.update(validated_data)
                except Exception as val_err:
                    emit('upload_error', {
                        'file_id': file_id,
                        'message': f"Validation Error: {getattr(val_err, 'message', str(val_err))}"
                    })
                    return
                
                state = self.upload_manager.initialize_upload(
                    file_id, filename, total_size, "project_creation", content_type, extra_data
                )
                
                if state:
                    ack_type = 'upload_start_ack'
                    message = 'Project creation upload started'
                    if state.get('status') == 'resuming':
                        ack_type = 'upload_resume_ack'
                        message = f"Resuming project creation upload from byte {state.get('received_bytes', 0)}"
                    
                    emit(ack_type, {
                        'file_id': file_id,
                        'received_bytes': state.get('received_bytes', 0),
                        'message': message
                    })
                    upload_log(f"{message} for file {filename} ({file_id})")
                else:
                    emit('upload_error', {'file_id': file_id, 'message': 'Failed to initialize project creation upload'})
                    
            except Exception as e:
                upload_log(f"Error in project creation upload_start: {e}")
                emit('upload_error', {'message': f'Upload start error: {str(e)}'})

        @self.socketio.on('upload_chunk', namespace=self.namespace)
        def handle_upload_chunk(data):
            """Handle file chunk upload for project creation."""
            try:
                file_id = data.get('file_id')
                chunk_num = int(data.get('chunk_num', 0))
                chunk_data_b64 = data.get('data')
                
                if not all([file_id, chunk_data_b64 is not None]):
                    emit('upload_error', {'message': 'Missing required chunk data'})
                    return
                
                try:
                    chunk_data = base64.b64decode(chunk_data_b64)
                except Exception:
                    emit('upload_error', {'file_id': file_id, 'message': 'Invalid chunk data encoding'})
                    return
                
                updated_state = self.upload_manager.write_chunk(file_id, chunk_data)
                if updated_state:
                    # Send progress updates periodically
                    if chunk_num % 5 == 0 or updated_state['received_bytes'] == updated_state['total_size']:
                        progress = min(100, int((updated_state['received_bytes'] / updated_state['total_size']) * 100)) if updated_state['total_size'] > 0 else 100
                        emit('upload_progress', {
                            'file_id': file_id,
                            'received': updated_state['received_bytes'],
                            'total': updated_state['total_size'],
                            'progress': progress
                        })
                else:
                    emit('upload_error', {'file_id': file_id, 'message': 'Failed to write chunk'})
                    
            except Exception as e:
                upload_log(f"Error in project creation upload_chunk: {e}")
                emit('upload_error', {'message': f'Chunk upload error: {str(e)}'})

        @self.socketio.on('upload_end', namespace=self.namespace)
        def handle_upload_end(data):
            """Handle upload completion for project creation."""
            try:
                file_id = data.get('file_id')
                if not file_id:
                    emit('upload_error', {'message': 'file_id required for upload_end'})
                    return
                
                final_state = self.upload_manager.finalize_upload(file_id)
                if final_state:
                    self.upload_manager.queue_file_for_processing(final_state, None, self.namespace)
                    
                    # Calculate upload stats
                    try:
                        start_time = datetime.fromisoformat(final_state.get('start_time_iso', datetime.now(timezone.utc).isoformat()))
                        elapsed = (datetime.now(timezone.utc) - start_time).total_seconds()
                        speed = final_state['total_size'] / (elapsed or 0.001) / 1024
                    except Exception:
                        elapsed = -1
                        speed = -1
                    
                    emit('upload_end_ack', {
                        'file_id': file_id,
                        'message': 'Project creation upload completed and verified. Queued for processing.',
                        'upload_speed_KBps': round(speed, 2) if speed != -1 else 'N/A',
                        'upload_time_sec': round(elapsed, 2) if elapsed != -1 else 'N/A'
                    })
                    upload_log(f"Project creation upload completed and queued: {final_state['filename']} ({file_id})")
                else:
                    emit('upload_error', {'file_id': file_id, 'message': 'Upload finalization failed'})
                    
            except Exception as e:
                upload_log(f"Error in project creation upload_end: {e}")
                emit('upload_error', {'message': f'Upload end error: {str(e)}'})

        @self.socketio.on('upload_cancel', namespace=self.namespace)
        def handle_upload_cancel(data):
            """Handle upload cancellation for project creation."""
            try:
                file_id = data.get('file_id')
                if not file_id:
                    emit('upload_error', {'message': 'file_id required for cancel'})
                    return
                
                upload_log(f"Cancelling project creation upload: {file_id}")
                self.upload_manager.cancel_upload(file_id)
                emit('upload_cancel_ack', {'file_id': file_id, 'message': 'Project creation upload cancelled'})
                
            except Exception as e:
                upload_log(f"Error in project creation upload_cancel: {e}")
                emit('upload_error', {'message': f'Cancel error: {str(e)}'})

        @self.socketio.on('ping', namespace=self.namespace)
        def handle_ping(data):
            """Handle ping for connection testing."""
            emit('pong', {'message': 'pong'})

class UploadSocketAPI:
    """SocketIO-based file upload API for a specific project."""
    
    def __init__(self, project: MDVProject, socketio: SocketIO, upload_manager: FileUploadManager):
        self.project = project
        self.socketio = socketio
        self.upload_manager = upload_manager
        self.namespace = f"/project/{project.id}"
        
        # Register event handlers
        self._register_upload_events()
        upload_log(f"UploadSocketAPI initialized for project {project.id}")

    def _register_upload_events(self):
        """Register SocketIO event handlers for upload operations."""
        
        @self.socketio.on('upload_query', namespace=self.namespace)
        def handle_upload_query(data):
            """Handle upload status query."""
            file_id = data.get('file_id')
            if not file_id:
                emit('upload_error', {'message': 'file_id required for query'})
                return
                
            upload_log(f"Upload query for file_id: {file_id}")
            state = self.upload_manager._load_upload_state(file_id)
            
            if state:
                status = state.get('status')
                if status in ['completed', 'queued', 'processing']:
                    upload_log(f"File {file_id} status: {status}, re-initiating processing")
                    self.upload_manager.queue_file_for_processing(state, None, self.namespace)
                    emit('upload_processing_initiated', {
                        'file_id': file_id,
                        'message': f'Processing for {state.get("original_filename", "file")} initiated'
                    })
                elif status in ['uploading', 'resuming', 'error']:
                    emit('upload_resume_info', {
                        'file_id': file_id,
                        'received_bytes': state.get('received_bytes', 0),
                        'total_size': state.get('total_size', 0)
                    })
                else:
                    emit('upload_not_found', {'file_id': file_id})
            else:
                emit('upload_not_found', {'file_id': file_id})

        @self.socketio.on('upload_start', namespace=self.namespace)
        def handle_upload_start(data):
            """Handle upload start."""
            try:
                file_id = data.get('file_id') or str(uuid.uuid4())
                filename = secure_filename(data.get('filename', f'file_{file_id}'))
                total_size = int(data.get('size', 0))
                content_type = data.get('content_type', 'application/octet-stream')
                
                extra_data = {}
                
                # Validation based on content type
                if content_type == "text/csv":
                    try:
                        validation_params = {k: data.get(k) for k in ["name", "view", "replace", "supplied_only"]}
                        validation_params = {k: v for k, v in validation_params.items() if v is not None}
                        with self.upload_manager.app.app_context():
                            validated_data = validate_datasource(self.project, validation_params)
                        extra_data.update(validated_data)
                    except Exception as val_err:
                        emit('upload_error', {
                            'file_id': file_id, 
                            'message': f"Validation Error: {getattr(val_err, 'message', str(val_err))}"
                        })
                        return
                elif content_type == "application/x-hdf" or filename.endswith('.h5ad'):
                    # Validate AnnData file
                    try:
                        validation_params = {}  # No specific params needed for AnnData
                        with self.upload_manager.app.app_context():
                            validated_data = validate_anndata(self.project, validation_params)
                        extra_data.update(validated_data)
                    except Exception as val_err:
                        emit('upload_error', {
                            'file_id': file_id,
                            'message': f"Validation Error: {getattr(val_err, 'message', str(val_err))}"
                        })
                        return
                
                state = self.upload_manager.initialize_upload(
                    file_id, filename, total_size, self.project.id, content_type, extra_data
                )
                
                if state:
                    ack_type = 'upload_start_ack'
                    message = 'Upload started'
                    if state.get('status') == 'resuming':
                        ack_type = 'upload_resume_ack'
                        message = f"Resuming upload from byte {state.get('received_bytes', 0)}"
                    
                    emit(ack_type, {
                        'file_id': file_id,
                        'received_bytes': state.get('received_bytes', 0),
                        'message': message
                    })
                    upload_log(f"{message} for file {filename} ({file_id})")
                else:
                    emit('upload_error', {'file_id': file_id, 'message': 'Failed to initialize upload'})
                    
            except Exception as e:
                upload_log(f"Error in upload_start: {e}")
                emit('upload_error', {'message': f'Upload start error: {str(e)}'})

        @self.socketio.on('upload_chunk', namespace=self.namespace)
        def handle_upload_chunk(data):
            """Handle file chunk upload."""
            try:
                file_id = data.get('file_id')
                chunk_num = int(data.get('chunk_num', 0))
                chunk_data_b64 = data.get('data')
                
                if not all([file_id, chunk_data_b64 is not None]):
                    emit('upload_error', {'message': 'Missing required chunk data'})
                    return
                
                try:
                    chunk_data = base64.b64decode(chunk_data_b64)
                except Exception:
                    emit('upload_error', {'file_id': file_id, 'message': 'Invalid chunk data encoding'})
                    return
                
                updated_state = self.upload_manager.write_chunk(file_id, chunk_data)
                if updated_state:
                    # Send progress updates periodically
                    if chunk_num % 5 == 0 or updated_state['received_bytes'] == updated_state['total_size']:
                        progress = min(100, int((updated_state['received_bytes'] / updated_state['total_size']) * 100)) if updated_state['total_size'] > 0 else 100
                        emit('upload_progress', {
                            'file_id': file_id,
                            'received': updated_state['received_bytes'],
                            'total': updated_state['total_size'],
                            'progress': progress
                        })
                else:
                    emit('upload_error', {'file_id': file_id, 'message': 'Failed to write chunk'})
                    
            except Exception as e:
                upload_log(f"Error in upload_chunk: {e}")
                emit('upload_error', {'message': f'Chunk upload error: {str(e)}'})

        @self.socketio.on('upload_end', namespace=self.namespace)
        def handle_upload_end(data):
            """Handle upload completion."""
            try:
                file_id = data.get('file_id')
                if not file_id:
                    emit('upload_error', {'message': 'file_id required for upload_end'})
                    return
                
                final_state = self.upload_manager.finalize_upload(file_id)
                if final_state:
                    self.upload_manager.queue_file_for_processing(final_state, None, self.namespace)
                    
                    # Calculate upload stats
                    try:
                        start_time = datetime.fromisoformat(final_state.get('start_time_iso', datetime.now(timezone.utc).isoformat()))
                        elapsed = (datetime.now(timezone.utc) - start_time).total_seconds()
                        speed = final_state['total_size'] / (elapsed or 0.001) / 1024
                    except Exception:
                        elapsed = -1
                        speed = -1
                    
                    emit('upload_end_ack', {
                        'file_id': file_id,
                        'message': 'Upload completed and verified. Queued for processing.',
                        'upload_speed_KBps': round(speed, 2) if speed != -1 else 'N/A',
                        'upload_time_sec': round(elapsed, 2) if elapsed != -1 else 'N/A'
                    })
                    upload_log(f"Upload completed and queued: {final_state['filename']} ({file_id})")
                else:
                    emit('upload_error', {'file_id': file_id, 'message': 'Upload finalization failed'})
                    
            except Exception as e:
                upload_log(f"Error in upload_end: {e}")
                emit('upload_error', {'message': f'Upload end error: {str(e)}'})

        @self.socketio.on('upload_cancel', namespace=self.namespace)
        def handle_upload_cancel(data):
            """Handle upload cancellation."""
            try:
                file_id = data.get('file_id')
                if not file_id:
                    emit('upload_error', {'message': 'file_id required for cancel'})
                    return
                
                upload_log(f"Cancelling upload: {file_id}")
                self.upload_manager.cancel_upload(file_id)
                emit('upload_cancel_ack', {'file_id': file_id, 'message': 'Upload cancelled'})
                
            except Exception as e:
                upload_log(f"Error in upload_cancel: {e}")
                emit('upload_error', {'message': f'Cancel error: {str(e)}'})

        @self.socketio.on('ping', namespace=self.namespace)
        def handle_ping(data):
            """Handle ping for connection testing."""
            emit('pong', {'message': 'pong'})

# --- Public API Functions ---

def initialize_socketio_upload(app: Flask, base_temp_dir: Optional[str] = None) -> 'FileUploadManager':
    """
    Initialize the SocketIO upload system. 
    This function should be called after mdv_socketio(app) has been called.
    
    Args:
        app (Flask): Flask application instance  
        base_temp_dir (Optional[str]): Base directory for temporary uploads
        
    Returns:
        FileUploadManager: Initialized upload manager
        
    Example:
        app = Flask(__name__)
        mdv_socketio(app)  # Initialize SocketIO first
        upload_manager = initialize_socketio_upload(app)
    """
    global _upload_manager, _socketio_instance, _project_creation_api
    
    if _upload_manager is not None:
        upload_log("Upload system already initialized")
        return _upload_manager
    
    # Get the SocketIO instance from the app extensions
    # This assumes mdv_socketio(app) was called first
    if hasattr(app, 'extensions') and 'socketio' in app.extensions:
        _socketio_instance = app.extensions['socketio']
    else:
        # Fallback: try to import from mdvtools.websocket
        try:
            from mdvtools.websocket import socketio as imported_socketio
            _socketio_instance = imported_socketio
            upload_log("Retrieved SocketIO instance from mdvtools.websocket")
        except ImportError:
            raise RuntimeError("SocketIO not found. Make sure mdv_socketio(app) is called before initialize_socketio_upload(app)")
    
    if _socketio_instance is None:
        raise RuntimeError("SocketIO instance not available. Call mdv_socketio(app) first.")
    
    _upload_manager = FileUploadManager(app, _socketio_instance, base_temp_dir)
    
    # Initialize project creation upload API for the main namespace
    _project_creation_api = ProjectCreationUploadAPI(app, _socketio_instance, _upload_manager)
    
    upload_log("SocketIO upload system initialized successfully")
    return _upload_manager

def register_project_for_upload(project_id: str, project: MDVProject) -> UploadSocketAPI:
    """
    Register a project for upload capabilities.
    
    Args:
        project_id (str): Unique project identifier
        project (MDVProject): Project instance
        
    Returns:
        UploadSocketAPI: Upload API instance for the project
        
    Example:
        register_project_for_upload(project.id, project)
    """
    global _upload_projects_map, _upload_manager, _socketio_instance
    
    if _upload_manager is None or _socketio_instance is None:
        raise RuntimeError("Upload system not initialized. Call initialize_socketio_upload(app) first.")
    
    _upload_projects_map[project_id] = project
    upload_api = UploadSocketAPI(project, _socketio_instance, _upload_manager)
    upload_log(f"Registered project {project_id} for uploads")
    return upload_api

def initialize_project_creation_uploads(app: Flask) -> bool:
    """
    Initialize project creation upload capability.
    This should be called during multi-project server setup.
    
    Args:
        app (Flask): Flask application instance with projects_base_dir configured
        
    Returns:
        bool: True if successfully initialized
        
    Example:
        app = Flask(__name__)
        app.config['projects_base_dir'] = '/path/to/projects'
        mdv_socketio(app)
        initialize_socketio_upload(app)
        initialize_project_creation_uploads(app)
    """
    global _project_creation_api
    
    if not app.config.get('projects_base_dir'):
        upload_log("Warning: projects_base_dir not configured, project creation uploads disabled")
        return False
        
    if _upload_manager is None or _socketio_instance is None:
        upload_log("Upload system not initialized, call initialize_socketio_upload() first")
        return False
        
    if _project_creation_api is not None:
        upload_log("Project creation uploads already initialized")
        return True
    
    # Initialize project creation upload API
    _project_creation_api = ProjectCreationUploadAPI(app, _socketio_instance, _upload_manager)
    upload_log("Project creation upload API initialized successfully")
    return True

# --- Legacy/Alternative API Functions ---

def initialize_upload_system(app: Flask, socketio: SocketIO, base_temp_dir: Optional[str] = None) -> FileUploadManager:
    """
    Legacy function name for backward compatibility.
    Use initialize_socketio_upload() instead.
    """
    global _upload_manager, _socketio_instance
    _socketio_instance = socketio
    _upload_manager = FileUploadManager(app, socketio, base_temp_dir)
    upload_log("Upload system initialized via legacy function")
    return _upload_manager

def get_upload_manager() -> Optional[FileUploadManager]:
    """Get the current upload manager instance."""
    return _upload_manager

def get_registered_projects() -> List[str]:
    """Get list of registered project IDs."""
    return list(_upload_projects_map.keys())

def remove_project_from_uploads(project_id: str) -> bool:
    """
    Remove a project from upload capabilities.
    
    Args:
        project_id (str): Project ID to remove
        
    Returns:
        bool: True if removed, False if not found
    """
    global _upload_projects_map
    if project_id in _upload_projects_map:
        del _upload_projects_map[project_id]
        upload_log(f"Removed project {project_id} from uploads")
        return True
    else:
        upload_log(f"Project {project_id} not found in upload registry")
        return False

def is_upload_system_initialized() -> bool:
    """Check if the upload system is initialized."""
    return _upload_manager is not None

def get_upload_stats() -> Dict[str, Any]:
    """Get upload system statistics."""
    if _upload_manager is None:
        return {"status": "not_initialized"}
    
    return {
        "status": "initialized",
        "base_temp_dir": _upload_manager.base_temp_dir,
        "registered_projects": get_registered_projects(),
        "upload_ttl_hours": _upload_manager.upload_ttl.total_seconds() / 3600,
        "processing_queue_size": len(_upload_manager.processing_queue) if _upload_manager else 0,
        "project_creation_enabled": _project_creation_api is not None
    }
