"""
Build information module for gathering version metadata from environment or git.

This module provides utilities to extract build/version information that can be
used for reproducibility and debugging. It checks environment variables first
(as set by Docker builds or devcontainers), then falls back to git commands.
"""
import os
import subprocess
from typing import Optional, Dict, Any
from datetime import datetime


def _get_from_env() -> Optional[Dict[str, Any]]:
    """Try to get build info from environment variables."""
    git_commit_hash = os.environ.get("GIT_COMMIT_HASH")
    git_commit_date = os.environ.get("GIT_COMMIT_DATE")
    git_branch = os.environ.get("GIT_BRANCH_NAME")
    build_date = os.environ.get("BUILD_DATE")
    git_dirty = os.environ.get("GIT_DIRTY", "false").lower() == "true"
    
    if git_commit_hash:
        return {
            "git_commit_hash": git_commit_hash,
            "git_commit_date": git_commit_date,
            "git_branch": git_branch,
            "build_date": build_date,
            "git_dirty": git_dirty,
            "source": "environment"
        }
    return None


def _get_from_git() -> Optional[Dict[str, Any]]:
    """Try to get build info from git commands."""
    try:
        # Start from the directory where this module is located
        git_root = os.path.dirname(__file__)
        
        # Verify we're in a git repository
        subprocess.run(
            ["git", "rev-parse", "--git-dir"],
            capture_output=True,
            check=True,
            timeout=5,
            cwd=git_root
        )
        
        # Get commit hash
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            cwd=git_root
        )
        if result.returncode != 0:
            return None
        git_commit_hash = result.stdout.strip()
        
        # Get commit date
        result = subprocess.run(
            ["git", "log", "-1", "--format=%cI"],
            capture_output=True,
            text=True,
            timeout=5,
            cwd=git_root
        )
        git_commit_date = result.stdout.strip() if result.returncode == 0 else None
        
        # Get branch name
        result = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            cwd=git_root
        )
        git_branch = result.stdout.strip() if result.returncode == 0 else None
        
        # Check if working directory is dirty
        result = subprocess.run(
            ["git", "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=5,
            cwd=git_root
        )
        git_dirty = len(result.stdout.strip()) > 0
        
        # Use current time as build date
        build_date = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
        
        return {
            "git_commit_hash": git_commit_hash,
            "git_commit_date": git_commit_date,
            "git_branch": git_branch,
            "build_date": build_date,
            "git_dirty": git_dirty,
            "source": "git"
        }
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError, OSError):
        return None


def get_build_info() -> Dict[str, Any]:
    """
    Get build information from environment variables or git.
    
    Returns:
        Dictionary with build information:
        - git_commit_hash: Full commit hash
        - git_commit_date: ISO format commit date
        - git_branch: Branch name
        - build_date: Build timestamp
        - git_dirty: Whether working directory has uncommitted changes
        - source: "environment", "git", or "unknown"
    """
    # Try environment first (Docker/build-time)
    info = _get_from_env()
    if info:
        return info
    
    # Fallback to git
    info = _get_from_git()
    if info:
        return info
    
    # No information available
    return {
        "git_commit_hash": None,
        "git_commit_date": None,
        "git_branch": None,
        "build_date": None,
        "git_dirty": None,
        "source": "unknown"
    }


def build_info_to_markdown(info: Dict[str, Any]) -> str:
    """
    Convert build info dictionary to markdown format.
    
    Args:
        info: Build info dictionary from get_build_info()
        
    Returns:
        Markdown formatted string
    """
    if info["source"] == "unknown":
        return "### Build information:\n\n*Build information not available.*\n"
    
    lines = ["### Build information:\n"]
    
    if info["git_commit_hash"]:
        lines.append(f"- **Commit hash**: `{info['git_commit_hash']}`")
    if info["git_commit_date"]:
        lines.append(f"- **Commit date**: {info['git_commit_date']}")
    if info["git_branch"]:
        lines.append(f"- **Branch**: {info['git_branch']}")
    if info["build_date"]:
        lines.append(f"- **Build date**: {info['build_date']}")
    if info["git_dirty"] is not None:
        lines.append(f"- **Working directory**: {'dirty' if info['git_dirty'] else 'clean'}")
    lines.append(f"- **Source**: {info['source']}")
    
    return "\n".join(lines) + "\n"


def get_build_info_markdown() -> str:
    """
    Convenience function to get build info and return it as markdown.
    
    Returns:
        Markdown formatted string with build information
    """
    return build_info_to_markdown(get_build_info())

