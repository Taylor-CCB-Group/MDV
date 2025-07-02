#!/bin/bash

# Safely extract git information without shell injection vulnerabilities
# This script should be used instead of direct git commands in shell contexts

set -euo pipefail

# Function to safely extract commit message
safe_commit_message() {
    # Use printf to avoid shell interpretation
    git show -s --format=%s | tr -d '\n\r' | sed 's/[^[:print:]]//g' | head -c 1000
}

# Function to safely extract commit hash
safe_commit_hash() {
    git rev-parse HEAD
}

# Function to safely extract branch name
safe_branch_name() {
    git rev-parse --abbrev-ref HEAD
}

# Function to safely extract commit date
safe_commit_date() {
    git log -1 --format=%cI
}

# Function to check if working directory is dirty
safe_git_dirty() {
    if [ -z "$(git status --porcelain)" ]; then
        echo "false"
    else
        echo "true"
    fi
}

# Main execution
case "${1:-}" in
    "commit_message")
        safe_commit_message
        ;;
    "commit_hash")
        safe_commit_hash
        ;;
    "branch_name")
        safe_branch_name
        ;;
    "commit_date")
        safe_commit_date
        ;;
    "git_dirty")
        safe_git_dirty
        ;;
    *)
        echo "Usage: $0 {commit_message|commit_hash|branch_name|commit_date|git_dirty}"
        exit 1
        ;;
esac 