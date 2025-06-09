from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional, Dict, Any
import json
from pathlib import Path
import os

@dataclass
class ChatLogItem:
    """Represents a single chat log entry"""
    context: str
    query: str
    prompt_template: str
    response: str
    timestamp: str
    conversation_id: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return {
            "context": self.context,
            "query": self.query,
            "prompt_template": self.prompt_template,
            "response": self.response,
            "timestamp": self.timestamp,
            "conversation_id": self.conversation_id
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ChatLogItem':
        """Create from dictionary (e.g. from JSON)"""
        return cls(
            context=data["context"],
            query=data["query"],
            prompt_template=data["prompt_template"],
            response=data["response"],
            timestamp=data["timestamp"],
            conversation_id=data.get("conversation_id")
        )

class ChatLogger:
    """Handles chat logging functionality"""
    def __init__(self, log_file_path: str):
        self.log_file_path = Path(log_file_path)
        self._ensure_log_file_exists()

    def _ensure_log_file_exists(self):
        """Ensure the log file exists and is properly initialized"""
        if not self.log_file_path.exists():
            self.log_file_path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.log_file_path, 'w') as f:
                json.dump([], f)

    def log_chat(self, item: ChatLogItem):
        """Log a chat item to the JSON file"""
        try:
            # Read existing logs
            with open(self.log_file_path, 'r') as f:
                logs = json.load(f)
            
            # Append new log
            logs.append(item.to_dict())
            
            # Write back
            with open(self.log_file_path, 'w') as f:
                json.dump(logs, f, indent=4)
        except Exception as e:
            print(f"Error logging chat: {e}")

    def get_logs(self) -> List[ChatLogItem]:
        """Get all chat logs"""
        try:
            with open(self.log_file_path, 'r') as f:
                logs = json.load(f)
            return [ChatLogItem.from_dict(log) for log in logs]
        except Exception as e:
            print(f"Error reading chat logs: {e}")
            return []

    def get_conversation_logs(self, conversation_id: str) -> List[ChatLogItem]:
        """Get logs for a specific conversation"""
        return [log for log in self.get_logs() if log.conversation_id == conversation_id] 