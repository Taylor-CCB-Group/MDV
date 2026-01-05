"""Simple MCP-Bio HTTP client for ChatMDV integration.

This module provides a client for communicating with the MCP-Bio server,
which provides bioinformatics analysis tools (UMAP, clustering, etc.).
"""
import httpx
from typing import Dict, Any, List, Optional
import os
import logging

logger = logging.getLogger(__name__)


class MCPBioClient:
    """HTTP client for MCP-Bio server communication, supporting standard MCP-over-SSE."""
    
    def __init__(self, server_url: str = "http://localhost:8000"):
        """Initialize the MCP-Bio client.
        
        Args:
            server_url: Base URL of the MCP-Bio server
        """
        self.server_url = server_url.rstrip("/")
        self._tools_cache: Optional[List[Dict]] = None
        self._client = httpx.Client(timeout=httpx.Timeout(300.0, connect=30.0))
    
    def health_check(self) -> bool:
        """Check if MCP-Bio server is reachable.
        
        Returns:
            True if server is healthy, False otherwise
        """
        try:
            # Try /health first
            response = self._client.get(f"{self.server_url}/health", timeout=2.0)
            if response.status_code == 200:
                return True
                
            # If 404, try /sse using stream to avoid hanging on the streaming response
            try:
                with self._client.stream("GET", f"{self.server_url}/sse", timeout=2.0) as response:
                    if response.status_code == 200:
                        return True
            except Exception:
                pass
                
            # Any response from uvicorn counts as "reachable"
            response = self._client.get(self.server_url, timeout=2.0)
            if "uvicorn" in response.headers.get("server", "").lower():
                return True
                
            return False
        except Exception:
            return False

    def _mcp_request(self, method: str, params: Dict[str, Any]) -> Dict[str, Any]:
        """Perform a standard MCP JSON-RPC request over SSE, including initialization."""
        request_id = f"req_{os.getpid()}_{method.replace('/', '_')}"
        
        try:
            # Use a fresh client for the SSE connection to avoid state issues
            with httpx.Client(timeout=httpx.Timeout(300.0, connect=30.0)) as client:
                # 1. Start SSE connection
                with client.stream("GET", f"{self.server_url}/sse") as sse_response:
                    sse_response.raise_for_status()
                    
                    msg_endpoint = None
                    lines_gen = sse_response.iter_lines()
                    
                    # 2. Read stream until we get the endpoint URL
                    for line in lines_gen:
                        if line.startswith("data: "):
                            path = line[len("data: "):].strip()
                            if path.startswith("/"):
                                msg_endpoint = f"{self.server_url}{path}"
                            else:
                                msg_endpoint = path
                            break
                    
                    if not msg_endpoint:
                        raise ConnectionError("Could not discover MCP message endpoint")
                    
                    # 3. Initialize handshake (required by standard MCP servers)
                    init_payload = {
                        "jsonrpc": "2.0",
                        "id": "init_auto",
                        "method": "initialize",
                        "params": {
                            "protocolVersion": "2024-11-05",
                            "capabilities": {},
                            "clientInfo": {"name": "ChatMDV", "version": "1.0.0"}
                        }
                    }
                    client.post(msg_endpoint, json=init_payload).raise_for_status()
                    
                    # Wait for any response (usually the init response)
                    for line in lines_gen:
                        if line.startswith("data: "):
                            break
                    
                    # 4. Send the actual POST request
                    payload = {
                        "jsonrpc": "2.0",
                        "id": request_id,
                        "method": method,
                        "params": params
                    }
                    post_response = client.post(msg_endpoint, json=payload)
                    post_response.raise_for_status()
                    
                    # 5. Continue reading the SSE stream until we find our response
                    for line in lines_gen:
                        if line.startswith("data: "):
                            data_str = line[len("data: "):].strip()
                            try:
                                import json
                                msg = json.loads(data_str)
                                # Check if this is the response to our request
                                if msg.get("id") == request_id:
                                    return msg
                            except Exception:
                                continue
                            
            raise ConnectionError(f"Stream closed before receiving response for {method}")
        except Exception as e:
            logger.error(f"MCP request {method} failed: {e}")
            raise

    def list_tools(self) -> List[Dict]:
        """List available tools from MCP-Bio server.
        
        Returns:
            List of tool definitions with name, description, and parameters
        """
        if self._tools_cache:
            return self._tools_cache
        
        try:
            # Try MCP-over-SSE first
            result = self._mcp_request("tools/list", {})
            
            if "error" in result:
                error_msg = result["error"].get("message", "Unknown error")
                logger.error(f"MCP list_tools error: {error_msg}")
                raise ConnectionError(f"MCP error: {error_msg}")
                
            tools = result.get("result", {}).get("tools", [])
            self._tools_cache = tools
            return self._tools_cache
        except Exception as e:
            logger.warning(f"Standard MCP list_tools failed, trying fallback: {e}")
            # Fallback for simple HTTP servers
            try:
                response = self._client.get(f"{self.server_url}/tools")
                response.raise_for_status()
                self._tools_cache = response.json()
                return self._tools_cache
            except Exception:
                logger.error(f"MCP list_tools fallback also failed")
                raise ConnectionError(f"Could not list tools from {self.server_url}") from e
    
    def call_tool(self, tool_name: str, arguments: Dict[str, Any]) -> Dict:
        """Call an MCP-Bio tool.
        
        Args:
            tool_name: Name of the tool to call (e.g., "umap", "leiden")
            arguments: Tool arguments as a dictionary
            
        Returns:
            Tool result as a dictionary with status and output data
        """
        try:
            # Try MCP-over-SSE first
            result = self._mcp_request("tools/call", {
                "name": tool_name,
                "arguments": arguments
            })
            
            if "error" in result:
                error_msg = result["error"].get("message", "Unknown error")
                logger.error(f"MCP call_tool {tool_name} error: {error_msg}")
                return {"status": "error", "error": error_msg}
                
            # Extract content from result
            content_list = result.get("result", {}).get("content", [])
            output_text = ""
            for item in content_list:
                if item.get("type") == "text":
                    output_text += item.get("text", "")
            
            # Attempt to parse as JSON if it looks like one
            try:
                import json
                parsed_output = json.loads(output_text)
                if isinstance(parsed_output, dict):
                    return parsed_output
            except Exception:
                pass
                
            return {
                "status": "success",
                "output": output_text
            }
            
        except Exception as e:
            logger.warning(f"Standard MCP call_tool failed, trying fallback: {e}")
            # Fallback for simple HTTP servers
            try:
                response = self._client.post(f"{self.server_url}/tools/{tool_name}", json=arguments)
                if response.status_code == 200:
                    return response.json()
            except Exception:
                pass
                    
            return {
                "status": "error", 
                "error": f"Tool call failed: {str(e)}"
            }
    
    def clear_cache(self):
        """Clear the tools cache to force refresh on next list_tools call."""
        self._tools_cache = None
    
    def __del__(self):
        """Clean up HTTP client on deletion."""
        if hasattr(self, '_client'):
            self._client.close()


# Singleton pattern for global client access
_client: Optional[MCPBioClient] = None


def get_mcp_client() -> MCPBioClient:
    """Get or create the global MCP-Bio client instance.
    
    The server URL is read from MCP_BIO_URL environment variable,
    defaulting to http://localhost:8000.
    
    Returns:
        MCPBioClient instance
    """
    global _client
    if _client is None:
        url = os.getenv("MCP_BIO_URL", "http://localhost:8000")
        
        # If we are in a container and using localhost, it might fail.
        # We can try to detect if host.docker.internal is available.
        if url == "http://localhost:8000":
            try:
                import socket
                # Try to resolve host.docker.internal
                socket.gethostbyname("host.docker.internal")
                url = "http://host.docker.internal:8000"
                logger.info(f"Detected host.docker.internal, using {url} as default")
            except socket.gaierror:
                # Fallback to localhost if host.docker.internal is not resolvable
                pass
                
        _client = MCPBioClient(url)
        logger.info(f"Initialized MCP-Bio client with URL: {url}")
    return _client


def reset_mcp_client():
    """Reset the global client (useful for testing or reconfiguration)."""
    global _client
    if _client is not None:
        _client._client.close()
    _client = None
