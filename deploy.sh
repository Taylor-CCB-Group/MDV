#!/bin/bash

# MDV Cross-Platform Setup Script
# Compatible with Windows WSL, macOS, and Linux

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logo
echo -e "${BLUE}"
echo "  __  __  _____   __      __"
echo " |  \/  | |  __ \ \ \    / /"
echo " | \  / | | |  | | \ \  / / "
echo " | |\/| | | |  | |  \ \/ /  "
echo " | |  | | | |__| |   \  /   "
echo " |_|  |_| |_____/     \/    "
echo ""
echo "Multi-Dimensional Viewer Setup"
echo -e "${NC}"

# Check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}
# Docker checks
check_docker() {
    echo -e "${BLUE}Checking Docker installation...${NC}"
    
    if ! command_exists docker; then
        echo -e "${RED}Error: Docker is not installed.${NC}"
        exit 1
    fi
    
    if ! docker info >/dev/null 2>&1; then
        echo -e "${RED}Error: Docker daemon is not running.${NC}"
        exit 1
    fi
    
    if ! command_exists docker-compose && ! docker compose version >/dev/null 2>&1; then
        echo -e "${RED}Error: Docker Compose is not available.${NC}"
        exit 1
    fi
    
    echo -e "${GREEN}✓ Docker is properly installed and running${NC}"
}

# Function to check for existing Docker containers, volumes, or networks
check_existing_service() {
  echo -e "${BLUE}Checking for existing MDV containers, volumes, or networks...${NC}"

  # Look for MDV containers, volumes, or networks
  existing_containers=$(docker ps -aq --filter "name=mdv")
  existing_volumes=$(docker volume ls -q --filter "name=mdv")
  existing_networks=$(docker network ls -q --filter "name=mdv")

  if [[ -n "$existing_containers" || -n "$existing_volumes" || -n "$existing_networks" ]]; then
    echo -e "${RED}✗ Error: Existing MDV Docker containers, volumes, or networks detected.${NC}"
    echo -e "${YELLOW}To remove them safely, run the following command:${NC}"
    echo
    echo -e "${BLUE}docker compose -p mdv down --remove-orphans${NC}"
    echo -e "${YELLOW}⚠ If you also want to delete volumes (this will erase database/data), run:${NC}"
    echo -e "${BLUE}docker compose -p mdv down -v --remove-orphans${NC}"
    echo
    exit 1
  else
    echo -e "${GREEN}✓ No conflicting MDV Docker containers, volumes, or networks found${NC}"
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

    # Ensure each prompt starts on a new line
    echo

    if [ "$hide_value" == "true" ]; then
        # Hide input (silent mode)
        read -s -p "Enter value for $var_name [$default_value]: " new_value
        echo "" >&2 # Move to the next line after input
    else
        # Show input with default value suggestion
        read -p "Enter value for $var_name [$default_value]: " new_value
        echo "" # Move to the next line after input
    fi

    new_value=${new_value:-$default_value}  # Use the new value or fallback to the current value
    echo "$var_name=$new_value"  # Return the value assignment (not the value alone)
  }

  # Check if Postgres volume exists
  DB_VOLUME_EXISTS=$(docker volume ls -q | grep -E ".*_postgres-data" | wc -l)

  if [ "$DB_VOLUME_EXISTS" -gt 0 ]; then
      echo -e "${YELLOW}⚠ Existing Postgres data volume detected."
      echo -e "DB_USER, DB_PASSWORD, and DB_NAME will NOT be changed to keep your database data safe.${NC}"

      # Preserve previously sourced values or use defaults if missing
      DB_USER="${DB_USER:-testuser}"
      DB_PASSWORD="${DB_PASSWORD:-testpass}"
      DB_NAME="${DB_NAME:-testdb}"
  else
      echo -e "${YELLOW}Enter DB credentials for new database setup:${NC}"
      echo -e "${YELLOW}If you want to keep DB data in future redeploys, keep these values the same.${NC}"

      # Capture the output of prompt_variable and extract the values
      DB_USER=$(prompt_variable "DB_USER" false "testuser"| cut -d'=' -f2 | tr -d '\n' | xargs)
      DB_PASSWORD=$(prompt_variable "DB_PASSWORD" true "testpass"| cut -d'=' -f2 | tr -d '\n' | xargs)
      DB_NAME=$(prompt_variable "DB_NAME" false "testdb"| cut -d'=' -f2 | tr -d '\n' | xargs)
  fi

  # Set PostgreSQL variables
  POSTGRES_USER=${DB_USER}
  POSTGRES_PASSWORD=${DB_PASSWORD}
  POSTGRES_DB=${DB_NAME}

  # Write new .env file
  {
    echo "# Flask Configuration"
    echo "FLASK_ENV=production"
    echo "PYTHONUNBUFFERED=1"

    echo "# Database Configuration"
    echo "DB_USER=$DB_USER"
    echo "DB_PASSWORD=$DB_PASSWORD"
    echo "DB_NAME=$DB_NAME"
    echo "DB_HOST=mdv_db"  # Default without prompting

    echo "POSTGRES_USER=$POSTGRES_USER"
    echo "POSTGRES_PASSWORD=$POSTGRES_PASSWORD"
    echo "POSTGRES_DB=$POSTGRES_DB"

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
      echo "$(prompt_variable "AUTH0_DB_CONNECTION" false)"
      echo "$(prompt_variable "AUTH0_PUBLIC_KEY_URI" false)"
      echo "$(prompt_variable "SHIBBOLETH_LOGIN_URL" false)"
      echo "$(prompt_variable "SHIBBOLETH_LOGOUT_URL" false)"
    else
      echo "ENABLE_AUTH=0"
    fi
    # Chat
    read -p "Enable Chat? (y/n): " enable_chat
    if [[ "$enable_chat" == "y" ]]; then
      echo "ENABLE_CHAT=1"
      echo "$(prompt_variable "OPENAI_API_KEY" false)"
    else
      echo "ENABLE_CHAT=0"
    fi
  } > "$env_file"

  echo -e "${GREEN}✓ Environment variables updated and saved to $env_file."
}

# Function to download and run Docker Compose
run_docker_compose() {
  local compose_url=$1
  local compose_file=$(basename "$compose_url")

  if [ ! -f "$compose_file" ]; then
    echo "Downloading $compose_file..."
    if ! curl -O "$compose_url"; then
      echo -e "${RED} Error: Failed to download $compose_file."
      exit 1
    fi
  fi

  echo "Pulling latest Docker images..."
  docker compose -f "$compose_file" pull

  echo "Starting Docker Compose..."
  docker compose -f "$compose_file" up -d || {
    echo -e "${RED} Error: Failed to start Docker Compose."
    exit 1
  }

  echo -e "${GREEN}✓ Docker Compose started successfully!"
}

# Function to open the application in the browser
open_browser() {
  URL="http://localhost:5055"
  echo -e "${GREEN}✓ Application deployed successfully!"
  echo -e "${GREEN} Open $URL in your browser."
}

# Main Execution
echo "Starting MDV Deployment..."

check_docker
#check_existing_service

# Create app_logs directory if it doesn't exist
#mkdir -p app_logs && echo "Ensured app_logs directory exists."

create_or_validate_env_file

DOCKER_COMPOSE_URL="https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/main/docker-local.yml"
run_docker_compose "$DOCKER_COMPOSE_URL"

open_browser
