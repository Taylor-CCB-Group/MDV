#!/bin/bash

# Function to check if Docker is installed
check_docker_installed() {
  if ! command -v docker &> /dev/null; then
    echo "Error: Docker is not installed. Please install Docker before running this script."
    exit 1
  fi
}

# Function to check if Docker daemon is running
check_docker_daemon() {
  if ! docker info &> /dev/null; then
    echo "Error: Docker daemon is not running. Please start Docker before running this script."
    exit 1
  fi
}

# Function to create or validate the .env file
create_or_validate_env_file() {
  env_file=".env"

  if [ ! -f "$env_file" ]; then
    echo "$env_file does not exist. Creating it now..."
    touch "$env_file"
  fi

  # Source existing .env file to retain values
  [ -f "$env_file" ] && source "$env_file"

  # Function to prompt user for input
  prompt_variable() {
    local var_name=$1
    local hide_value=$2
    local default_value=$3
    local current_value=${!var_name:-$default_value}

    read -p "Enter value for $var_name [$current_value]: " new_value
    new_value=${new_value:-$current_value}

    echo "$var_name=$new_value"
  }

  # Write new .env file
  {
    echo "# Flask Configuration"
    echo "FLASK_ENV=production"
    echo "PYTHONUNBUFFERED=1"

    echo "# Database Configuration"
    echo "$(prompt_variable "DB_USER" false)"
    echo "$(prompt_variable "DB_PASSWORD" true)"
    echo "$(prompt_variable "DB_NAME" false)"
    echo "DB_HOST=mdv_db"  # Default without prompting

    echo "POSTGRES_USER=$DB_USER"
    echo "POSTGRES_PASSWORD=$DB_PASSWORD"
    echo "POSTGRES_DB=$DB_NAME"

    # Authentication
    read -p "Enable Authentication? (y/n): " enable_auth
    if [[ "$enable_auth" == "y" ]]; then
      echo "ENABLE_AUTH=1"
      echo "$(prompt_variable "FLASK_SECRET_KEY" true)"
      echo "$(prompt_variable "LOGIN_REDIRECT_URL" false)"
      echo "$(prompt_variable "AUTH0_DOMAIN" false)"
      echo "$(prompt_variable "AUTH0_CLIENT_ID" false)"
      echo "$(prompt_variable "AUTH0_CLIENT_SECRET" true)"
      echo "$(prompt_variable "AUTH0_CALLBACK_URL" false)"
      echo "$(prompt_variable "AUTH0_AUDIENCE" false)"
      echo "$(prompt_variable "AUTH0_PUBLIC_KEY_URI" false)"
      echo "$(prompt_variable "SHIBBOLETH_LOGIN_URL" false)"
      echo "$(prompt_variable "SHIBBOLETH_LOGOUT_URL" false)"
    else
      echo "ENABLE_AUTH=0"
    fi
  } > "$env_file"

  echo "Environment variables updated and saved to $env_file."
}

# Function to download and run Docker Compose
run_docker_compose() {
  local compose_url=$1
  local compose_file=$(basename "$compose_url")

  if [ ! -f "$compose_file" ]; then
    echo "Downloading $compose_file..."
    if ! curl -O "$compose_url"; then
      echo "Error: Failed to download $compose_file."
      exit 1
    fi
  fi

  echo "Pulling latest Docker images..."
  docker compose -f "$compose_file" pull

  echo "Starting Docker Compose..."
  docker compose -f "$compose_file" up -d || {
    echo "Error: Failed to start Docker Compose."
    exit 1
  }

  echo "Docker Compose started successfully!"
}

# Function to open the application in the browser
open_browser() {
  URL="https://localhost:5055"
  echo "Application deployed successfully!"
  echo "Open $URL in your browser."
}

# Main Execution
echo "Starting MDV Deployment..."

check_docker_installed
check_docker_daemon

# Create app_logs directory if it doesn't exist
mkdir -p app_logs && echo "Ensured app_logs directory exists."

create_or_validate_env_file

DOCKER_COMPOSE_URL="https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/auth_rbac/docker-secrets-local.yml"
run_docker_compose "$DOCKER_COMPOSE_URL"

open_browser
