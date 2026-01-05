"""LangChain tool wrappers for MCP-Bio bioinformatics tools.

This module provides LangChain-compatible tool wrappers for the MCP-Bio server,
allowing the chat agent to run bioinformatics analyses like UMAP and clustering.
"""
import logging
from typing import Optional, List

from langchain.tools import tool
from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

# Module-level state for project context
_project = None
_h5ad_manager = None


def init_mcp_tools(project):
    """Initialize MCP tools with project reference.
    
    Must be called before using any MCP tools.
    
    Args:
        project: MDVProject instance
    """
    global _project, _h5ad_manager
    from mdvtools.mdvproject import MDVProject
    
    if not isinstance(project, MDVProject):
        raise TypeError(f"Expected MDVProject, got {type(project)}")
    
    _project = project
    
    # Lazy init of h5ad manager
    from .h5ad_manager import H5ADManager
    _h5ad_manager = H5ADManager(project)
    
    logger.info(f"Initialized MCP tools for project: {project.dir}")


def _check_initialized():
    """Check that MCP tools have been initialized."""
    if _project is None or _h5ad_manager is None:
        raise RuntimeError(
            "MCP tools not initialized. Call init_mcp_tools(project) first."
        )


def _get_path_variations(full_path: str) -> list:
    """Generate path variations to try for MCP-Bio file access.
    
    This handles the case where ChatMDV and MCP-Bio have different volume mounts.
    """
    import os
    
    filename = os.path.basename(full_path)
    rel_path = full_path.replace("/app/mdv/", "")
    
    host_mdv_path = os.getenv("MCP_HOST_MDV_PATH", "")
    
    path_variations = []
    
    # Priority 1: Translated host path (most likely to work)
    if host_mdv_path:
        host_full_path = os.path.join(host_mdv_path, rel_path)
        path_variations.append(host_full_path)
    
    # Priority 2: Try common macOS/Linux home expansions
    for home_base in ["/Users/stephent/mdv", "/home/pn/mdv", os.path.expanduser("~/mdv")]:
        candidate = os.path.join(home_base, rel_path)
        if candidate not in path_variations:
            path_variations.append(candidate)
    
    # Priority 3: Container paths (if MCP is also in Docker with shared volume)
    path_variations.extend([
        full_path, 
        rel_path, 
        f"./{rel_path}", 
        f"/app/shared/{rel_path}",
        f"/app/shared/{filename}",
        f"/data/{rel_path}",
        f"/shared/{rel_path}"
    ])
    
    return path_variations


def _is_path_error(result: dict) -> tuple[bool, str]:
    """Check if result indicates a path/file not found error.
    
    Returns (is_error, error_message).
    """
    error_msg = None
    if "error" in result:
        error_msg = str(result["error"])
    elif "output" in result and isinstance(result["output"], str) and "Error" in result["output"]:
        error_msg = result["output"]
    
    if error_msg:
        if "File not found" in error_msg or "No such file" in error_msg or "Unable to" in error_msg:
            return True, error_msg
    return False, ""


@tool("list_analysis_tools")
def list_analysis_tools() -> str:
    """List available bioinformatics analysis tools from MCP-Bio server.
    
    Use this tool to see what analyses can be run on the data.
    """
    from .mcp_client import get_mcp_client
    
    try:
        client = get_mcp_client()
        tools = client.list_tools()
        # #region agent log
        log_debug("H3", "Tools listed from MCP server", {"tools_count": len(tools), "tool_names": [t.get('name') for t in tools]})
        # #endregion
        
        if not tools:
            return "No analysis tools available. Is the MCP-Bio server running?"
        
        result = "**Available Analysis Tools:**\n\n"
        for t in tools:
            name = t.get('name', 'unknown')
            desc = t.get('description', 'No description')
            result += f"• **{name}**: {desc}\n"
        
        return result
        
    except ConnectionError as e:
        return f"**Error**: Could not connect to MCP-Bio server. {str(e)}"
    except Exception as e:
        logger.error(f"Error listing tools: {e}")
        return f"**Error**: {str(e)}"


class UMAPParams(BaseModel):
    """Parameters for UMAP dimensionality reduction."""
    n_neighbors: int = Field(
        default=15, 
        description="Number of neighbors for UMAP. Higher values preserve more global structure."
    )
    min_dist: float = Field(
        default=0.5, 
        description="Minimum distance between points. Lower values create tighter clusters."
    )


@tool("run_umap", args_schema=UMAPParams)
def run_umap(n_neighbors: int = 15, min_dist: float = 0.5) -> str:
    """Run UMAP dimensionality reduction on the cells.
    
    Use this tool when the user asks for:
    - UMAP
    - Dimensionality reduction
    - 2D embedding
    - Visualization of cell relationships
    
    After running, tells the user to reload the page to see results.
    """
    _check_initialized()
    from .mcp_client import get_mcp_client
    
    try:
        client = get_mcp_client()
        base_path = _h5ad_manager.get_path()
        path_variations = _get_path_variations(base_path)
        
        logger.info(f"Running UMAP: n_neighbors={n_neighbors}, min_dist={min_dist}")
        
        # #region agent log
        log_debug("H11", "run_umap path variations", {"paths": path_variations})
        # #endregion
        
        last_error = None
        for attempt_path in path_variations:
            output_path = attempt_path.replace(".h5ad", "_umap_result.h5ad")
            
            # #region agent log
            log_debug("H11", f"run_umap attempting", {"path": attempt_path})
            # #endregion
            
            try:
                result = client.call_tool("umap", {
                    "input_path": attempt_path,
                    "output_path": output_path,
                    "n_neighbors": n_neighbors,
                    "min_dist": min_dist
                })
            except Exception as call_err:
                # #region agent log
                log_debug("H11", f"run_umap exception", {"path": attempt_path, "error": str(call_err)})
                # #endregion
                last_error = str(call_err)
                continue  # Try next path
            
            is_path_err, err_msg = _is_path_error(result)
            if is_path_err:
                # #region agent log
                log_debug("H11", f"run_umap path error", {"path": attempt_path, "error": err_msg})
                # #endregion
                last_error = err_msg
                continue  # Try next path
            
            if result.get("status") == "success" or "output" in result:
                # Sync results back to MDV
                sync_result = _h5ad_manager.sync_to_mdv(output_path, ["umap_0", "umap_1"])
                
                synced = sync_result.get("synced", [])
                if synced:
                    return (
                        "**UMAP Complete!** ✓\n\n"
                        f"Parameters: n_neighbors={n_neighbors}, min_dist={min_dist}\n\n"
                        f"Added columns: `{'`, `'.join(synced)}`\n\n"
                        "**→ Reload the page to see the new columns in your charts.**"
                    )
                else:
                    return (
                        "**UMAP completed** but columns could not be synced to MDV.\n"
                        f"Result file: {output_path}"
                    )
            else:
                error = result.get("error", "Unknown error")
                return f"**UMAP Failed:** {error}"
        
        # All paths failed
        return f"**UMAP Failed:** Could not find file. Last error: {last_error}"
            
    except ConnectionError as e:
        return f"**Error**: Could not connect to MCP-Bio server. {str(e)}"
    except Exception as e:
        logger.error(f"UMAP error: {e}")
        return f"**Error running UMAP:** {str(e)}"


class ClusterParams(BaseModel):
    """Parameters for Leiden clustering."""
    resolution: float = Field(
        default=1.0, 
        description="Clustering resolution. Higher values create more clusters."
    )


@tool("run_clustering", args_schema=ClusterParams)
def run_clustering(resolution: float = 1.0) -> str:
    """Run Leiden clustering to group similar cells.
    
    Use this tool when the user asks for:
    - Clustering
    - Grouping cells
    - Finding cell populations
    - Leiden clustering
    - Cell type identification
    
    After running, tells the user to reload the page to see results.
    """
    _check_initialized()
    from .mcp_client import get_mcp_client
    
    try:
        client = get_mcp_client()
        base_path = _h5ad_manager.get_path()
        path_variations = _get_path_variations(base_path)
        
        logger.info(f"Running clustering: resolution={resolution}")
        
        last_error = None
        for attempt_path in path_variations:
            output_path = attempt_path.replace(".h5ad", "_clustered_result.h5ad")
            
            try:
                result = client.call_tool("leiden", {
                    "input_path": attempt_path,
                    "output_path": output_path,
                    "resolution": resolution
                })
            except Exception as call_err:
                last_error = str(call_err)
                continue  # Try next path
            
            is_path_err, err_msg = _is_path_error(result)
            if is_path_err:
                last_error = err_msg
                continue  # Try next path
            
            if result.get("status") == "success" or "n_clusters" in result:
                # Sync results back to MDV
                sync_result = _h5ad_manager.sync_to_mdv(output_path, ["leiden"])
                
                n_clusters = result.get("n_clusters", "?")
                synced = sync_result.get("synced", [])
                
                if synced:
                    return (
                        "**Clustering Complete!** ✓\n\n"
                        f"Found **{n_clusters} clusters** (resolution={resolution})\n\n"
                        f"Added column: `{'`, `'.join(synced)}`\n\n"
                        "**→ Reload the page to see the new column in your charts.**"
                    )
                else:
                    return (
                        f"**Clustering completed** with {n_clusters} clusters, "
                        "but the column could not be synced to MDV.\n"
                        f"Result file: {output_path}"
                    )
            else:
                error = result.get("error", "Unknown error")
                return f"**Clustering Failed:** {error}"
        
        # All paths failed
        return f"**Clustering Failed:** Could not find file. Last error: {last_error}"
            
    except ConnectionError as e:
        return f"**Error**: Could not connect to MCP-Bio server. {str(e)}"
    except Exception as e:
        logger.error(f"Clustering error: {e}")
        return f"**Error running clustering:** {str(e)}"


def get_mcp_tools() -> List:
    """Get all MCP tools for inclusion in the agent.
    
    Returns:
        List of LangChain tool functions
    """
    from .mcp_client import get_mcp_client
    from pydantic import create_model, Field as PydanticField
    from typing import Optional, Any
    import json

    tools = [list_analysis_tools, run_umap, run_clustering]
    
    try:
        client = get_mcp_client()
        mcp_tools = client.list_tools()
        
        # We already have these explicitly wrapped
        existing_names = ["umap", "leiden", "list_analysis_tools"]
        
        for t in mcp_tools:
            name = t.get('name')
            if name in existing_names:
                continue
                
            # Create a dynamic LangChain tool for each MCP tool
            desc = t.get('description', 'No description')
            schema = t.get('inputSchema', {})
            
            # Use a closure to capture the tool name
            def make_tool_func(tool_name, tool_schema, tool_desc):
                # Extract properties from schema and create a dynamic Pydantic model
                properties = tool_schema.get("properties", {})
                required = tool_schema.get("required", [])
                
                # Build field definitions for Pydantic create_model
                field_definitions = {}
                for prop_name, prop_info in properties.items():
                    prop_type = prop_info.get("type", "string")
                    prop_desc = prop_info.get("description", "")
                    
                    # Map JSON schema types to Python types
                    type_map = {
                        "string": str,
                        "integer": int,
                        "number": float,
                        "boolean": bool,
                        "array": list,
                        "object": dict
                    }
                    python_type = type_map.get(prop_type, str)
                    
                    # Make optional if not required
                    if prop_name in required:
                        field_definitions[prop_name] = (python_type, PydanticField(description=prop_desc))
                    else:
                        field_definitions[prop_name] = (Optional[python_type], PydanticField(default=None, description=prop_desc))
                
                # Create dynamic Pydantic model for this tool's parameters
                if field_definitions:
                    DynamicSchema = create_model(f"{tool_name}_params", **field_definitions)
                else:
                    DynamicSchema = None
                
                @tool(tool_name, args_schema=DynamicSchema)
                def dynamic_tool(**kwargs) -> str:
                    """Dynamic tool wrapper for MCP."""
                    _check_initialized()
                    
                    try:
                        full_path = _h5ad_manager.get_path()
                        # Variations to try:
                        # 1. Host path translation (if MCP_HOST_MDV_PATH is set)
                        # 2. Absolute path (/app/mdv/2/file.h5ad) - Best for shared volumes
                        # 3. Relative path from mdv root (2/file.h5ad)
                        # 4. Prefixed relative path (./2/file.h5ad)
                        # 5. Shared path as per docs (/app/shared/2/file.h5ad)
                        # 6. Data path (/data/2/file.h5ad)
                        
                        import os
                        filename = os.path.basename(full_path)
                        rel_path = full_path.replace("/app/mdv/", "")
                        
                        # Check if host path translation is configured
                        # This is needed when MCP-Bio runs on the host (not in Docker)
                        # and ChatMDV runs in a container with /app/mdv mounted from host's ~/mdv
                        host_mdv_path = os.getenv("MCP_HOST_MDV_PATH", "")
                        
                        path_variations = []
                        
                        # Priority 1: Translated host path (most likely to work)
                        if host_mdv_path:
                            host_full_path = os.path.join(host_mdv_path, rel_path)
                            path_variations.append(host_full_path)
                        
                        # Priority 2: Try common macOS/Linux home expansions
                        # These are guesses based on typical docker-compose setups
                        for home_base in ["/Users/stephent/mdv", "/home/pn/mdv", os.path.expanduser("~/mdv")]:
                            candidate = os.path.join(home_base, rel_path)
                            if candidate not in path_variations:
                                path_variations.append(candidate)
                        
                        # Priority 3: Container paths (if MCP is also in Docker with shared volume)
                        path_variations.extend([
                            full_path, 
                            rel_path, 
                            f"./{rel_path}", 
                            f"/app/shared/{rel_path}",
                            f"/app/shared/{filename}",
                            f"/data/{rel_path}",
                            f"/shared/{rel_path}"
                        ])
                        
                        # #region agent log
                        log_debug("H10", f"Dynamic tool {tool_name} preparing to try paths", {"variations": path_variations})
                        # #endregion

                        last_result = None
                        
                        for attempt_path in path_variations:
                            # Prepare kwargs for this attempt
                            current_kwargs = kwargs.copy()
                            
                            # Inject path into expected arguments if not provided by agent
                            if "file_path" in properties and ("file_path" not in current_kwargs or not current_kwargs["file_path"]):
                                current_kwargs["file_path"] = attempt_path
                            if "input_path" in properties and ("input_path" not in current_kwargs or not current_kwargs["input_path"]):
                                current_kwargs["input_path"] = attempt_path
                            if "output_path" in properties and ("output_path" not in current_kwargs or not current_kwargs["output_path"]):
                                current_kwargs["output_path"] = attempt_path.replace(".h5ad", f"_{tool_name}_result.h5ad")

                            # FIXED: Pass ALL kwargs through - don't filter out agent-provided parameters.
                            # The MCP server will validate. Only path injection happens above.
                            filtered_kwargs = current_kwargs

                            # #region agent log
                            log_debug("H10", f"Attempting {tool_name}", {"path": attempt_path, "kwargs": filtered_kwargs})
                            # #endregion
                            
                            try:
                                result = client.call_tool(tool_name, filtered_kwargs)
                                last_result = result
                            except Exception as e:
                                # #region agent log
                                log_debug("H10", f"Tool {tool_name} CRASHED", {"path": attempt_path, "error": str(e)})
                                # #endregion
                                last_result = {"error": str(e)}
                                continue
                            
                            # Check if it succeeded
                            if isinstance(result, dict):
                                # MCP can return errors in two ways:
                                # 1. {"error": "..."} - explicit error key
                                # 2. {"status": "success", "output": "Error calling tool..."} - error buried in output
                                error_msg = None
                                if "error" in result:
                                    error_msg = str(result["error"])
                                elif "output" in result and isinstance(result["output"], str) and "Error" in result["output"]:
                                    # Check if output contains an error message (e.g. file not found)
                                    error_msg = result["output"]
                                
                                if error_msg is None:
                                    # #region agent log
                                    log_debug("H10", f"Tool {tool_name} SUCCESS", {"path": attempt_path, "result": result})
                                    # #endregion
                                    return json.dumps(result, indent=2)
                                else:
                                    # #region agent log
                                    log_debug("H10", f"Tool {tool_name} FAILED", {"path": attempt_path, "error": error_msg})
                                    # #endregion
                                    if "File not found" in error_msg or "No such file" in error_msg or "Unable to" in error_msg:
                                        # Keep trying other paths
                                        continue
                                    else:
                                        # It's a different kind of error (e.g. parameter error), return it
                                        return json.dumps(result, indent=2)
                            else:
                                # Not a dict, maybe just a string success or failure
                                # #region agent log
                                log_debug("H10", f"Tool {tool_name} returned non-dict", {"path": attempt_path, "result": str(result)})
                                # #endregion
                                return str(result)

                        # If we get here, all variations failed
                        return json.dumps(last_result or {"error": "All path variations failed"}, indent=2)
                        
                    except Exception as e:
                        return f"Error: {str(e)}"
                
                # Update the docstring to match the MCP tool's description
                dynamic_tool.description = f"{tool_name}: {tool_desc}"
                return dynamic_tool
            
            tools.append(make_tool_func(name, schema, desc))
            
    except Exception as e:
        logger.warning(f"Could not load dynamic MCP tools: {e}")
        
    return tools


def is_mcp_available() -> bool:
    """Check if MCP-Bio server is available.
    
    Returns:
        True if server is reachable, False otherwise
    """
    try:
        from .mcp_client import get_mcp_client
        client = get_mcp_client()
        return client.health_check()
    except Exception:
        return False


# #region agent log
import json
import time
def log_debug(hypothesis_id, message, data=None):
    log_entry = {
        "id": f"log_{int(time.time())}_{hypothesis_id}",
        "timestamp": int(time.time() * 1000),
        "location": "mcp_tools.py",
        "message": message,
        "data": data or {},
        "sessionId": "debug-session",
        "runId": "run1",
        "hypothesisId": hypothesis_id
    }
    try:
        with open("/app/.cursor/debug.log", "a") as f:
            f.write(json.dumps(log_entry) + "\n")
    except Exception:
        pass
# #endregion
