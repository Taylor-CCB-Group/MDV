import { io, type Socket } from 'socket.io-client';

export interface SocketIOUploadConfig {
    namespace: string;
    file: File;
    fileName?: string;
    fileId?: string;
    chunkSize?: number;
    maxRetries?: number;
    retryDelay?: number;
    // CSV-specific options
    datasourceName?: string;
    replace?: boolean;
    view?: string;
    suppliedOnly?: boolean;
    // Progress callback
    onProgress?: (progress: number, uploaded: number, total: number) => void;
    // Status callbacks
    onStatusChange?: (status: string, message?: string) => void;
    onError?: (error: any) => void;
    onSuccess?: (result: any) => void;
    socketPath?: string;
}

export interface SocketIOUploadState {
    fileId: string;
    status: 'initializing' | 'uploading' | 'processing' | 'completed' | 'error' | 'cancelled';
    uploadedBytes: number;
    totalBytes: number;
    error?: string;
    result?: any;
    socket?: Socket;
}

export class SocketIOUploadClient {
    private config: SocketIOUploadConfig;
    private state: SocketIOUploadState;
    private socket: Socket | null = null;
    private connectionEstablished = false;
    private uploadAcknowledged = false;
    private serverRespondedToQuery = false;
    private uploadTransferComplete = false;
    private processingComplete = false;
    private shouldExitProcessingWait = false;
    private resumeOffset = 0;
    private isResuming = false;
    private cancelRequested = false;

    constructor(config: SocketIOUploadConfig) {
        this.config = {
            chunkSize: 256 * 1024, // 256KB
            maxRetries: 3,
            retryDelay: 5000, // 5 seconds
            ...config
        };

        this.state = {
            fileId: config.fileId || this.generateFileId(),
            status: 'initializing',
            uploadedBytes: 0,
            totalBytes: config.file.size,
        };
    }

    private generateFileId(): string {
        return `upload_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`;
    }

    private detectFileType(): { type: string; contentType: string } {
        const fileName = this.config.fileName || this.config.file.name;
        const lowerName = fileName.toLowerCase();
        
        if (lowerName.endsWith('.csv')) {
            return { type: 'csv', contentType: 'text/csv' };
        } else if (lowerName.endsWith('.h5ad')) {
            return { type: 'anndata', contentType: 'application/x-hdf' };
        } else if (lowerName.endsWith('.tiff') || lowerName.endsWith('.tif')) {
            return { type: 'tiff', contentType: 'image/tiff' };
        } else if (lowerName.endsWith('.zip')) {
            return { type: 'zip', contentType: 'application/zip' };
        }
        
        return { type: 'unknown', contentType: 'application/octet-stream' };
    }

    private updateStatus(status: SocketIOUploadState['status'], message?: string) {
        this.state.status = status;
        if (this.config.onStatusChange) {
            this.config.onStatusChange(status, message);
        }
    }

    private updateProgress(uploaded: number, total: number) {
        this.state.uploadedBytes = uploaded;
        this.state.totalBytes = total;
        const progress = total > 0 ? (uploaded / total) * 100 : 0;
        
        if (this.config.onProgress) {
            this.config.onProgress(progress, uploaded, total);
        }
    }

    private setupSocketEventHandlers() {
        if (!this.socket) return;

        this.socket.on('connect', () => {
            console.log('Connected to SocketIO server');
            this.connectionEstablished = true;
        });

        this.socket.on('disconnect', () => {
            console.log('Disconnected from SocketIO server');
            this.connectionEstablished = false;
        });

        this.socket.on('upload_start_ack', (data: any) => {
            console.log('Upload start acknowledged:', data);
            if (data.file_id === this.state.fileId) {
                this.uploadAcknowledged = true;
            }
        });

        this.socket.on('upload_resume_ack', (data: any) => {
            console.log('Upload resume acknowledged:', data);
            if (data.file_id === this.state.fileId) {
                this.resumeOffset = data.received_bytes || 0;
                this.isResuming = true;
                this.uploadAcknowledged = true;
            }
        });

        this.socket.on('upload_end_ack', (data: any) => {
            console.log('Upload end acknowledged:', data);
            if (data.file_id === this.state.fileId) {
                this.uploadTransferComplete = true;
            }
        });

        this.socket.on('upload_progress', (data: any) => {
            if (data.file_id === this.state.fileId) {
                this.updateProgress(data.received || 0, data.total || this.state.totalBytes);
            }
        });

        this.socket.on('upload_processing_initiated', (data: any) => {
            console.log('Server initiated processing:', data);
            if (data.file_id === this.state.fileId) {
                this.updateStatus('processing', 'Server processing file...');
                this.serverRespondedToQuery = true;
            }
        });

        this.socket.on('upload_processing', (data: any) => {
            console.log('Server processing file:', data);
            if (data.file_id === this.state.fileId) {
                this.updateStatus('processing', 'File processing started');
            }
        });

        this.socket.on('upload_success', (data: any) => {
            console.log('Processing successful:', data);
            if (data.file_id === this.state.fileId) {
                this.processingComplete = true;
                this.shouldExitProcessingWait = true;
                this.state.result = data;
                this.updateStatus('completed', 'File processed successfully');
                if (this.config.onSuccess) {
                    this.config.onSuccess(data);
                }
            }
        });

        this.socket.on('upload_error', (data: any) => {
            console.error('Server error:', data);
            if (data.file_id === this.state.fileId || !data.file_id) {
                this.processingComplete = true;
                this.shouldExitProcessingWait = true;
                this.state.error = data.message || 'Unknown error';
                this.updateStatus('error', this.state.error);
                if (this.config.onError) {
                    this.config.onError(data);
                }
            }
        });

        this.socket.on('upload_resume_info', (data: any) => {
            console.log('Resume info received:', data);
            if (data.file_id === this.state.fileId) {
                this.resumeOffset = data.received_bytes || 0;
                this.isResuming = true;
                this.uploadTransferComplete = false;
                this.serverRespondedToQuery = true;
            }
        });

        this.socket.on('upload_not_found', (data: any) => {
            console.log('Server does not have state for file:', data);
            if (data.file_id === this.state.fileId) {
                this.resumeOffset = 0;
                this.isResuming = false;
                this.uploadTransferComplete = false;
                this.serverRespondedToQuery = true;
            }
        });

        this.socket.on('pong', () => {
            console.log('Received pong from server');
        });
    }

    private async connectToServer(): Promise<boolean> {
        return new Promise((resolve) => {
            console.log('Using namespace:', this.config.namespace);
            console.log('Using socket path:', this.config.socketPath);
            
            // Connect to the specific namespace directly
            this.socket = io(`${this.config.namespace}`, {
                path: this.config.socketPath || '/socket.io',   // '/test/socket.io' or '/carroll/socket.io'
                autoConnect: false,
                transports: ['polling'],
                timeout: 60000,
                forceNew: true,
                reconnection: false,
                reconnectionDelay: 1000,
                reconnectionAttempts: 5,
            });

            this.setupSocketEventHandlers();

            this.socket.connect();

            // Wait for connection with timeout
            const connectionTimeout = setTimeout(() => {
                if (!this.connectionEstablished) {
                    console.error('Connection timeout after 15 seconds');
                    resolve(false);
                }
            }, 15000);

            this.socket.on('connect', () => {
                console.log('Socket connected successfully to namespace:', this.config.namespace);
                clearTimeout(connectionTimeout);
                this.connectionEstablished = true;
                resolve(true);
            });

            this.socket.on('connect_error', (error) => {
                console.error('Socket connection error:', error);
                clearTimeout(connectionTimeout);
                resolve(false);
            });

            this.socket.on('disconnect', (reason) => {
                console.log('Socket disconnected:', reason);
                this.connectionEstablished = false;
            });
        });
    }


    private async queryUploadStatus(): Promise<boolean> {
        return new Promise((resolve) => {
            if (!this.socket) {
                resolve(false);
                return;
            }

            // Don't query if we're already in processing state
            if (this.state.status === 'processing' && this.processingComplete === false) {
                console.log('Already in processing state, skipping query');
                resolve(true);
                return;
            }

            this.serverRespondedToQuery = false;
            const queryMsg = { file_id: this.state.fileId };
            
            console.log('Querying server status for file_id:', this.state.fileId);
            this.socket.emit('upload_query', queryMsg);

            // Wait for response with timeout
            const timeout = setTimeout(() => {
                if (!this.serverRespondedToQuery) {
                    console.log('Query timeout, assuming new upload');
                    this.resumeOffset = 0;
                    this.isResuming = false;
                    resolve(true);
                }
            }, 20000);

            // Check for response in a polling manner
            const checkResponse = () => {
                if (this.serverRespondedToQuery) {
                    clearTimeout(timeout);
                    resolve(true);
                } else {
                    setTimeout(checkResponse, 500);
                }
            };
            checkResponse();
        });
    }

    private async startUpload(): Promise<void> {
        if (!this.socket) throw new Error('Socket not connected');

        const { type, contentType } = this.detectFileType();
        
        const startMsg: any = {
            file_id: this.state.fileId,
            filename: this.config.fileName || this.config.file.name,
            size: this.config.file.size,
            content_type: contentType,
        };

        // Add CSV-specific parameters
        if (type === 'csv') {
            startMsg.name = this.config.datasourceName || 
                           this.config.fileName?.replace(/\.[^/.]+$/, '') || 
                           this.config.file.name.replace(/\.[^/.]+$/, '');
            startMsg.view = this.config.view || 'default';
            startMsg.replace = this.config.replace || false;
            startMsg.supplied_only = this.config.suppliedOnly || false;
        }

        console.log('Sending upload start for', type, 'file, ID', this.state.fileId);
        this.socket.emit('upload_start', startMsg);
    }

    private async sendFileChunks(): Promise<void> {
        if (!this.socket) throw new Error('Socket not connected');

        const chunkSize = this.config.chunkSize || 256 * 1024;
        let chunkNum = 0;
        let uploadedBytes = this.resumeOffset;
        let currentOffset = this.resumeOffset;

        console.log('Starting file transmission from offset', this.resumeOffset);

        while (!this.cancelRequested && currentOffset < this.config.file.size) {
            const endOffset = Math.min(currentOffset + chunkSize, this.config.file.size);
            const chunk = this.config.file.slice(currentOffset, endOffset);

            try {
                // Read chunk as base64
                const base64Data = await this.readFileChunkAsBase64(chunk);
                
                const chunkMsg = {
                    file_id: this.state.fileId,
                    chunk_num: chunkNum,
                    data: base64Data,
                };

                // Send chunk
                this.socket.emit('upload_chunk', chunkMsg);
                
                uploadedBytes += chunk.size;
                currentOffset = endOffset;
                chunkNum++;

                // Update progress
                this.updateProgress(uploadedBytes, this.config.file.size);

                // Small delay to prevent overwhelming the server
                await this.delay(10);
                
            } catch (error) {
                console.error('Error processing chunk:', error);
                throw error;
            }
        }

        console.log('Finished sending chunks. Total bytes uploaded:', uploadedBytes);
    }

    private readFileChunkAsBase64(chunk: Blob): Promise<string> {
        return new Promise((resolve, reject) => {
            const reader = new FileReader();
            reader.onload = () => {
                const dataUrl = reader.result as string;
                resolve(dataUrl.substring(dataUrl.indexOf(',') + 1));
            };
            reader.onerror = (error) => {
                reject(error);
            };
            reader.readAsDataURL(chunk);
        });
    }
    
    private delay(ms: number): Promise<void> {
        return new Promise(resolve => setTimeout(resolve, ms));
    }

    private async endUpload(): Promise<void> {
        if (!this.socket) throw new Error('Socket not connected');

        const endMsg = { file_id: this.state.fileId };
        console.log('Sending upload end for file_id:', this.state.fileId);
        this.socket.emit('upload_end', endMsg);
    }

private async waitForProcessing(): Promise<void> {
    const timeout = 3600000; // 1 hour timeout
    const startTime = Date.now();
    let lastPingTime = Date.now();
    const pingInterval = 30000; // 30 seconds
    let reconnectAttempts = 0;
    const maxReconnectAttempts = 5;

    console.log('Waiting for processing completion...');

    while (!this.shouldExitProcessingWait && !this.cancelRequested) {
        const currentTime = Date.now();

        // Check timeout
        if (currentTime - startTime > timeout) {
            throw new Error('Timeout waiting for processing completion');
        }

        // Check connection and handle reconnection
        if (!this.socket?.connected) {
            console.warn('Connection lost during processing wait');
            
            if (reconnectAttempts < maxReconnectAttempts) {
                try {
                    console.log(`Attempting reconnection ${reconnectAttempts + 1}/${maxReconnectAttempts}`);
                    
                    // Disconnect old socket completely
                    if (this.socket) {
                        this.socket.disconnect();
                        this.socket = null;
                    }
                    
                    // Create new connection with preferred transport
                    const connected = await this.connectToServer();
                    if (connected) {
                        console.log('Reconnected successfully');
                        reconnectAttempts = 0; // Reset on successful connection
                        lastPingTime = currentTime; // Reset ping timer
                        
                        // Re-query status to get current processing state
                        await this.queryUploadStatus();
                    } else {
                        reconnectAttempts++;
                        console.warn(`Reconnection attempt ${reconnectAttempts} failed`);
                        
                        // Exponential backoff
                        await new Promise(resolve => 
                            setTimeout(resolve, Math.min(1000 * (2 ** reconnectAttempts), 30000))
                        );
                    }
                } catch (error) {
                    reconnectAttempts++;
                    console.warn('Reconnection failed:', error);
                    
                    if (reconnectAttempts >= maxReconnectAttempts) {
                        throw new Error(`Failed to reconnect after ${maxReconnectAttempts} attempts`);
                    }
                }
            } else {
                throw new Error(`Connection lost and failed to reconnect after ${maxReconnectAttempts} attempts`);
            }
        }

        // Send periodic pings only when connected
        if (this.socket?.connected && currentTime - lastPingTime > pingInterval) {
            try {
                this.socket.emit('ping', { message: 'keepalive', file_id: this.state.fileId });
                lastPingTime = currentTime;
                console.log('Sent keepalive ping');
            } catch (error) {
                console.warn('Failed to send keepalive ping:', error);
            }
        }

        // Wait before next check
        await new Promise(resolve => setTimeout(resolve, 1000));
    }
}

    public async upload(): Promise<any> {
        try {
            this.updateStatus('initializing', 'Connecting to server...');
            console.log('Starting upload for file:', this.config.file.name, 'Size:', this.config.file.size);

            // Connect to server
            const connected = await this.connectToServer();
            if (!connected) {
                throw new Error('Failed to connect to server');
            }

            // Query upload status
            const shouldProceed = await this.queryUploadStatus();
            if (!shouldProceed) {
                throw new Error('Server refused upload');
            }

            // Start upload if needed
            if (!this.uploadTransferComplete) {
                this.updateStatus('uploading', 'Starting upload...');
                
                await this.startUpload();
                console.log('Upload start sent, waiting for acknowledgment...');
                
                // Wait for acknowledgment
                await this.waitForAcknowledgment();
                console.log('Upload acknowledged, starting file transfer...');
                
                // Send file chunks
                await this.sendFileChunks();
                console.log('File chunks sent, ending upload...');
                
                // End upload
                await this.endUpload();
                console.log('Upload end sent, waiting for completion...');
                
                // Wait for upload completion acknowledgment
                await this.waitForUploadCompletion();
                console.log('Upload transfer completed');
            }

            // Wait for processing
            this.updateStatus('processing', 'Processing file...');
            console.log('Waiting for server processing...');
            await this.waitForProcessing();

            if (this.state.result) {
                console.log('Upload completed successfully:', this.state.result);
                return this.state.result;
            } else {
                throw new Error('Upload completed but no result received');
            }

        } catch (error) {
            console.error('Upload error:', error);
            this.state.error = error instanceof Error ? error.message : 'Unknown error';
            this.updateStatus('error', this.state.error);
            if (this.config.onError) {
                this.config.onError(error);
            }
            throw error;
        } finally {
            this.disconnect();
        }
    }

    private async waitForAcknowledgment(): Promise<void> {
        return new Promise((resolve, reject) => {
            const timeout = setTimeout(() => {
                reject(new Error('Upload acknowledgment timeout'));
            }, 15000);

            const checkAcknowledgment = () => {
                if (this.uploadAcknowledged) {
                    clearTimeout(timeout);
                    resolve();
                } else {
                    setTimeout(checkAcknowledgment, 500);
                }
            };
            checkAcknowledgment();
        });
    }

    private async waitForUploadCompletion(): Promise<void> {
        return new Promise((resolve, reject) => { // Add reject
            const timeout = setTimeout(() => {
                // FIX: Reject with an error
                reject(new Error('Timeout: Server did not acknowledge upload completion. The connection was likely lost.'));
            }, 30000);

            const checkCompletion = () => {
                if (this.cancelRequested) {
                    clearTimeout(timeout);
                    reject(new Error('Upload was cancelled.'));
                    return;
                }
                if (this.uploadTransferComplete) {
                    clearTimeout(timeout);
                    resolve();
                } else {
                    setTimeout(checkCompletion, 500);
                }
            };
            checkCompletion();
        });
    }

    public cancel(): void {
        this.cancelRequested = true;
        this.updateStatus('cancelled', 'Upload cancelled');
        
        if (this.socket?.connected) {
            this.socket.emit('upload_cancel', { file_id: this.state.fileId });
        }
        
        this.disconnect();
    }

    private disconnect(): void {
        if (this.socket) {
            this.socket.disconnect();
            this.socket = null;
        }
        this.connectionEstablished = false;
    }

    public getState(): SocketIOUploadState {
        return { ...this.state };
    }
}

// Factory function for easier usage
export function createSocketIOUpload(config: SocketIOUploadConfig): SocketIOUploadClient {
    return new SocketIOUploadClient(config);
}
