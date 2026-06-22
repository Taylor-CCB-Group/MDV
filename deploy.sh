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
DEPLOYMENT_NAME="mdv"
DEPLOY_MODE="replace"
APP_PORT="5055"
ENV_FILE=".env"
POSTGRES_MODE="dedicated"
REUSE_DB_CONTAINER="mdv-db"
REUSE_DB_NETWORK=""
DB_SCHEMA=""
DB_NAME="mdv"
SHARED_DB_NETWORK="mdv-shared-db"
SHARED_DB_HOST_PATH="./mdv-data/postgres-shared"
STORAGE_MODE="create"
HOST_DATA_PATH=""
RESET_DATA_VOLUME="0"
PREVIOUS_STORAGE_MODE=""
PREVIOUS_HOST_DATA_PATH=""
REMOVE_OLD_HOST_PATH="0"
REMOVE_OLD_VOLUME="0"
PREV_DB_BACKEND=""
PREV_DB_HOST=""
PREV_DB_SCHEMA=""
PREV_DB_USER=""
PREV_DB_NAME=""

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

prompt_deployment_settings() {
  get_next_available_port() {
    local port=5055
    while docker ps --format "{{.Ports}}" | grep -E "(0.0.0.0|\\[::\\]):${port}->" >/dev/null; do
      port=$((port + 1))
    done
    echo "$port"
  }
  get_existing_deployment_port() {
    local ports_line=""
    ports_line=$(docker ps \
      --filter "label=com.docker.compose.project=${DEPLOYMENT_NAME}" \
      --filter "label=com.docker.compose.service=mdv_app" \
      --format "{{.Ports}}" | head -n 1)
    local port_pattern='(0\.0\.0\.0|\[::\]|:::):([0-9]+)->5055'
    if [[ "$ports_line" =~ $port_pattern ]]; then
      echo "${BASH_REMATCH[2]}"
    fi
  }
  get_existing_storage_mount() {
    local app_container="${DEPLOYMENT_NAME}-mdv_app-1"
    local mount_type=""
    local mount_source=""
    if docker container inspect "$app_container" >/dev/null 2>&1; then
      mount_type=$(docker inspect -f '{{range .Mounts}}{{if eq .Destination "/app/mdv"}}{{.Type}}{{end}}{{end}}' "$app_container")
      mount_source=$(docker inspect -f '{{range .Mounts}}{{if eq .Destination "/app/mdv"}}{{.Source}}{{end}}{{end}}' "$app_container")
      if [[ -n "$mount_type" ]]; then
        echo "${mount_type}|${mount_source}"
      fi
    fi
  }

  read -p "Deployment name [mdv]: " deployment_name_input
  DEPLOYMENT_NAME=${deployment_name_input:-mdv}

  if docker ps -aq --filter "label=com.docker.compose.project=${DEPLOYMENT_NAME}" | grep -q .; then
    read -p "Deployment '${DEPLOYMENT_NAME}' exists. Action (redeploy/rename/exit) [redeploy]: " deploy_action_input
    DEPLOY_ACTION=$(echo "${deploy_action_input:-redeploy}" | tr '[:upper:]' '[:lower:]' | xargs)
    if [[ "$DEPLOY_ACTION" == "exit" ]]; then
      echo "Exiting deployment."
      exit 0
    elif [[ "$DEPLOY_ACTION" == "rename" ]]; then
      read -p "New deployment name: " new_deployment_name
      if [[ -z "$new_deployment_name" ]]; then
        echo -e "${RED}Error: deployment name cannot be empty.${NC}"
        exit 1
      fi
      DEPLOYMENT_NAME="$new_deployment_name"
      if docker ps -aq --filter "label=com.docker.compose.project=${DEPLOYMENT_NAME}" | grep -q .; then
        echo -e "${RED}Error: deployment '${DEPLOYMENT_NAME}' already exists. Choose a different name.${NC}"
        exit 1
      fi
      DEPLOY_MODE="new"
    else
      DEPLOY_MODE="replace"
    fi
  else
    DEPLOY_MODE="new"
    echo -e "${YELLOW}No existing deployment found for '${DEPLOYMENT_NAME}'. Using new mode.${NC}"
  fi

  DEFAULT_STORAGE_MODE="create"
  PREVIOUS_STORAGE_MOUNT=$(get_existing_storage_mount)
  if [[ -n "$PREVIOUS_STORAGE_MOUNT" ]]; then
    PREVIOUS_STORAGE_TYPE=${PREVIOUS_STORAGE_MOUNT%%|*}
    PREVIOUS_STORAGE_SOURCE=${PREVIOUS_STORAGE_MOUNT#*|}
    if [[ "$PREVIOUS_STORAGE_TYPE" == "bind" ]]; then
      PREVIOUS_STORAGE_MODE="host"
      PREVIOUS_HOST_DATA_PATH="$PREVIOUS_STORAGE_SOURCE"
      DEFAULT_STORAGE_MODE="host"
    elif [[ "$PREVIOUS_STORAGE_TYPE" == "volume" ]]; then
      PREVIOUS_STORAGE_MODE="reuse"
      DEFAULT_STORAGE_MODE="reuse"
    fi
  elif [[ "$DEPLOY_MODE" == "replace" ]] && docker volume inspect "${DEPLOYMENT_NAME}_mdv-data" >/dev/null 2>&1; then
    DEFAULT_STORAGE_MODE="reuse"
    PREVIOUS_STORAGE_MODE="reuse"
  fi

  APP_PORT_DEFAULT=$(get_next_available_port)
  if [[ "$DEPLOY_MODE" == "replace" ]]; then
    EXISTING_APP_PORT=$(get_existing_deployment_port)
    if [[ -n "$EXISTING_APP_PORT" ]]; then
      APP_PORT_DEFAULT="$EXISTING_APP_PORT"
      echo -e "${YELLOW}Detected existing app port ${APP_PORT_DEFAULT} for '${DEPLOYMENT_NAME}'.${NC}"
    fi
  fi
  read -p "Application host port [${APP_PORT_DEFAULT}]: " app_port_input
  APP_PORT=${app_port_input:-$APP_PORT_DEFAULT}

  ENV_FILE=".env.${DEPLOYMENT_NAME}"
}


# Function to create or validate the .env file
create_or_validate_env_file() {
  env_file="$ENV_FILE"

  if [ ! -f "$env_file" ]; then
    echo "$env_file does not exist. Creating it now..."
    touch "$env_file"
  fi

  # Source existing .env file to retain values
  [ -f "$env_file" ] && source "$env_file"
  PREV_DB_BACKEND="${DB_BACKEND:-}"
  PREV_DB_HOST="${DB_HOST:-}"
  PREV_DB_SCHEMA="${DB_SCHEMA:-}"
  PREV_DB_USER="${DB_USER:-}"
  PREV_DB_NAME="${DB_NAME:-}"

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

  ensure_shared_postgres() {
    if docker container inspect "$REUSE_DB_CONTAINER" >/dev/null 2>&1; then
      is_running=$(docker inspect -f '{{.State.Running}}' "$REUSE_DB_CONTAINER" 2>/dev/null || echo "false")
      if [[ "$is_running" != "true" ]]; then
        echo -e "${YELLOW}Shared postgres container '$REUSE_DB_CONTAINER' exists but is stopped. Starting it...${NC}"
        docker start "$REUSE_DB_CONTAINER" >/dev/null
      fi
      REUSE_DB_NETWORK=$(docker inspect "$REUSE_DB_CONTAINER" --format '{{range $k, $v := .NetworkSettings.Networks}}{{$k}}{{end}}')
      return 0
    fi

    read -p "Shared postgres container '$REUSE_DB_CONTAINER' not found. Create it now? (Y/n): " create_shared_input
    create_shared_choice=$(echo "${create_shared_input:-y}" | tr '[:upper:]' '[:lower:]' | xargs)
    if [[ "$create_shared_choice" == "n" || "$create_shared_choice" == "no" ]]; then
      echo -e "${RED}Cannot continue in reuse mode without a shared postgres container.${NC}"
      exit 1
    fi

    read -p "Shared Postgres host folder [${SHARED_DB_HOST_PATH}]: " shared_db_path_input
    SHARED_DB_HOST_PATH=${shared_db_path_input:-$SHARED_DB_HOST_PATH}
    docker network inspect "$SHARED_DB_NETWORK" >/dev/null 2>&1 || docker network create "$SHARED_DB_NETWORK" >/dev/null
    mkdir -p "$SHARED_DB_HOST_PATH"
    if ! docker run -d \
      --name "$REUSE_DB_CONTAINER" \
      --network "$SHARED_DB_NETWORK" \
      -e POSTGRES_USER="$DB_USER" \
      -e POSTGRES_PASSWORD="$DB_PASSWORD" \
      -e POSTGRES_DB="$DB_NAME" \
      -v "$SHARED_DB_HOST_PATH:/var/lib/postgresql/data" \
      postgres:16 >/dev/null; then
      echo -e "${RED}Failed to create shared postgres container '$REUSE_DB_CONTAINER'.${NC}"
      exit 1
    fi

    for _i in $(seq 1 30); do
      if docker exec "$REUSE_DB_CONTAINER" pg_isready -U "$DB_USER" -d "$DB_NAME" >/dev/null 2>&1; then
        REUSE_DB_NETWORK="$SHARED_DB_NETWORK"
        DB_HOST="$REUSE_DB_CONTAINER"
        return 0
      fi
      sleep 1
    done
    echo -e "${RED}Shared postgres container '$REUSE_DB_CONTAINER' did not become ready in time.${NC}"
    exit 1
  }

  # Select DB backend (default sqlite)
  echo
  read -p "Choose database backend (sqlite/postgres) [sqlite]: " db_backend_choice
  DB_BACKEND=${db_backend_choice:-${DB_BACKEND:-sqlite}}
  DB_BACKEND=$(echo "$DB_BACKEND" | tr '[:upper:]' '[:lower:]' | xargs)
  if [[ "$DB_BACKEND" != "sqlite" && "$DB_BACKEND" != "postgres" ]]; then
      echo -e "${YELLOW}Unknown backend '$DB_BACKEND'. Falling back to sqlite.${NC}"
      DB_BACKEND="sqlite"
  fi

  DB_HOST="mdv-db"
  SQLITE_DB_PATH="${SQLITE_DB_PATH:-/app/mdv/mdv.sqlite3}"
  DB_SCHEMA="$DEPLOYMENT_NAME"

  if [[ "$DB_BACKEND" == "postgres" ]]; then
    if docker container inspect "$REUSE_DB_CONTAINER" >/dev/null 2>&1; then
      REUSE_DB_NETWORK=$(docker inspect "$REUSE_DB_CONTAINER" --format '{{range $k, $v := .NetworkSettings.Networks}}{{$k}}{{end}}')
    fi
    EXISTING_POSTGRES_COUNT=0
    if docker container inspect "$REUSE_DB_CONTAINER" >/dev/null 2>&1; then
      EXISTING_POSTGRES_COUNT=1
    fi
    POSTGRES_MODE="reuse"
    echo -e "${YELLOW}Using shared postgres deployment mode (reuse).${NC}"
    DB_USER=$(prompt_variable "DB_USER" false "testuser" | cut -d'=' -f2 | tr -d '\n' | xargs)
    DB_PASSWORD=$(prompt_variable "DB_PASSWORD" true "testpass" | cut -d'=' -f2 | tr -d '\n' | xargs)
    DB_NAME=$(prompt_variable "DB_NAME" false "$DB_NAME" | cut -d'=' -f2 | tr -d '\n' | xargs)
    DB_HOST="$REUSE_DB_CONTAINER"
    ensure_shared_postgres
    DB_HOST=$(prompt_variable "DB_HOST" false "${DB_HOST:-$REUSE_DB_CONTAINER}" | cut -d'=' -f2 | tr -d '\n' | xargs)
  else
    # Keep defaults present in .env for compatibility if user switches backend later.
    DB_USER="${DB_USER:-testuser}"
    DB_PASSWORD="${DB_PASSWORD:-testpass}"
    DB_NAME="${DB_NAME:-mdv}"
    DB_SCHEMA=""
  fi

  if [[ "$DB_BACKEND" == "postgres" ]]; then
    STORAGE_MODE="reuse"
    echo -e "${YELLOW}App data storage prompt skipped for postgres backend.${NC}"
    if [[ "${PREV_DB_BACKEND:-}" == "sqlite" ]]; then
      if [[ "$PREVIOUS_STORAGE_MODE" == "host" && -n "$PREVIOUS_HOST_DATA_PATH" ]]; then
        read -p "Switching sqlite->postgres. Remove old host folder '${PREVIOUS_HOST_DATA_PATH}' and its data? (y/N): " remove_old_host_input
        remove_old_host_choice=$(echo "${remove_old_host_input:-n}" | tr '[:upper:]' '[:lower:]' | xargs)
        if [[ "$remove_old_host_choice" == "y" || "$remove_old_host_choice" == "yes" ]]; then
          REMOVE_OLD_HOST_PATH="1"
        fi
      fi
      if docker volume inspect "${DEPLOYMENT_NAME}_mdv-data" >/dev/null 2>&1; then
        read -p "Switching sqlite->postgres. Remove old app volume '${DEPLOYMENT_NAME}_mdv-data'? (y/N): " remove_old_volume_input
        remove_old_volume_choice=$(echo "${remove_old_volume_input:-n}" | tr '[:upper:]' '[:lower:]' | xargs)
        if [[ "$remove_old_volume_choice" == "y" || "$remove_old_volume_choice" == "yes" ]]; then
          REMOVE_OLD_VOLUME="1"
        fi
      fi
    fi
  else
    if [[ "$DEPLOY_MODE" == "replace" ]]; then
      if [[ "$DEFAULT_STORAGE_MODE" == "host" ]]; then
        echo -e "${YELLOW}Detected existing host-mapped data folder. Defaulting storage mode to host.${NC}"
      elif [[ "$DEFAULT_STORAGE_MODE" == "reuse" ]]; then
        echo -e "${YELLOW}Detected existing volume-based data mount. Defaulting storage mode to reuse.${NC}"
      fi
    fi
    read -p "Data storage mode (create/reuse/host) [${DEFAULT_STORAGE_MODE}]: " storage_mode_input
    STORAGE_MODE=$(echo "${storage_mode_input:-$DEFAULT_STORAGE_MODE}" | tr '[:upper:]' '[:lower:]' | xargs)
    if [[ "$STORAGE_MODE" != "create" && "$STORAGE_MODE" != "reuse" && "$STORAGE_MODE" != "host" ]]; then
      echo -e "${YELLOW}Unknown storage mode '$STORAGE_MODE'. Falling back to ${DEFAULT_STORAGE_MODE}.${NC}"
      STORAGE_MODE="$DEFAULT_STORAGE_MODE"
    fi
    if [[ "$STORAGE_MODE" == "host" ]]; then
      HOST_PATH_DEFAULT="./mdv-data/${DEPLOYMENT_NAME}"
      if [[ "$PREVIOUS_STORAGE_MODE" == "host" && -n "$PREVIOUS_HOST_DATA_PATH" ]]; then
        HOST_PATH_DEFAULT="$PREVIOUS_HOST_DATA_PATH"
      fi
      read -p "Host folder for /app/mdv [${HOST_PATH_DEFAULT}]: " host_data_input
      HOST_DATA_PATH=${host_data_input:-$HOST_PATH_DEFAULT}
      if [[ "$PREVIOUS_STORAGE_MODE" == "reuse" ]]; then
        read -p "Switching from named volume to host folder. Remove old volume '${DEPLOYMENT_NAME}_mdv-data'? (y/N): " remove_old_volume_input
        remove_old_volume_choice=$(echo "${remove_old_volume_input:-n}" | tr '[:upper:]' '[:lower:]' | xargs)
        if [[ "$remove_old_volume_choice" == "y" || "$remove_old_volume_choice" == "yes" ]]; then
          REMOVE_OLD_VOLUME="1"
        fi
      fi
    elif [[ "$STORAGE_MODE" == "create" && "$PREVIOUS_STORAGE_MODE" == "host" && -n "$PREVIOUS_HOST_DATA_PATH" ]]; then
      read -p "Switching from host folder to named volume. Remove old host folder '${PREVIOUS_HOST_DATA_PATH}' and its data? (y/N): " remove_old_host_input
      remove_old_host_choice=$(echo "${remove_old_host_input:-n}" | tr '[:upper:]' '[:lower:]' | xargs)
      if [[ "$remove_old_host_choice" == "y" || "$remove_old_host_choice" == "yes" ]]; then
        REMOVE_OLD_HOST_PATH="1"
      fi
    fi
    if [[ "$DEPLOY_MODE" == "replace" && "$STORAGE_MODE" == "create" ]] && docker volume inspect "${DEPLOYMENT_NAME}_mdv-data" >/dev/null 2>&1; then
      read -p "Existing data volume '${DEPLOYMENT_NAME}_mdv-data' found. Remove existing data before redeploy? (y/N): " remove_data_input
      remove_data_choice=$(echo "${remove_data_input:-n}" | tr '[:upper:]' '[:lower:]' | xargs)
      if [[ "$remove_data_choice" == "y" || "$remove_data_choice" == "yes" ]]; then
        RESET_DATA_VOLUME="1"
      else
        echo -e "${YELLOW}Keeping existing data volume (reusing existing data).${NC}"
        STORAGE_MODE="reuse"
      fi
    fi
  fi

  if [[ "$STORAGE_MODE" == "host" ]]; then
    if [[ "$PREVIOUS_STORAGE_MODE" == "host" && -n "$PREVIOUS_HOST_DATA_PATH" && "$PREVIOUS_HOST_DATA_PATH" != "$HOST_DATA_PATH" ]]; then
      read -p "Previous host folder '${PREVIOUS_HOST_DATA_PATH}' differs. Remove old folder and its data? (y/N): " remove_old_host_input
      remove_old_host_choice=$(echo "${remove_old_host_input:-n}" | tr '[:upper:]' '[:lower:]' | xargs)
      if [[ "$remove_old_host_choice" == "y" || "$remove_old_host_choice" == "yes" ]]; then
        REMOVE_OLD_HOST_PATH="1"
      fi
    fi
    if [[ -e "$HOST_DATA_PATH" && ! -d "$HOST_DATA_PATH" ]]; then
      echo -e "${RED}Error: '$HOST_DATA_PATH' exists but is not a folder.${NC}"
      exit 1
    fi
    if [[ -d "$HOST_DATA_PATH" ]] && [[ -n "$(ls -A "$HOST_DATA_PATH" 2>/dev/null)" ]]; then
      read -p "Host folder '$HOST_DATA_PATH' is not empty. Remove existing data before deployment? (y/N): " clear_host_input
      clear_host_choice=$(echo "${clear_host_input:-n}" | tr '[:upper:]' '[:lower:]' | xargs)
      if [[ "$clear_host_choice" == "y" || "$clear_host_choice" == "yes" ]]; then
        rm -rf "$HOST_DATA_PATH"/* "$HOST_DATA_PATH"/.[!.]* "$HOST_DATA_PATH"/..?* 2>/dev/null || true
      else
        echo -e "${YELLOW}Keeping existing host folder data.${NC}"
      fi
    fi
  fi

  POSTGRES_USER=${DB_USER}
  POSTGRES_PASSWORD=${DB_PASSWORD}
  POSTGRES_DB=${DB_NAME}

  # Write new .env file
  {
    echo "# Flask Configuration"
    echo "FLASK_ENV=production"
    echo "PYTHONUNBUFFERED=1"

    echo "# Database Configuration"
    echo "DB_BACKEND=$DB_BACKEND"
    echo "POSTGRES_MODE=$POSTGRES_MODE"
    echo "SQLITE_DB_PATH=$SQLITE_DB_PATH"
    echo "DB_USER=$DB_USER"
    echo "DB_PASSWORD=$DB_PASSWORD"
    echo "DB_NAME=$DB_NAME"
    echo "DB_HOST=$DB_HOST"
    echo "DB_SCHEMA=$DB_SCHEMA"

    echo "POSTGRES_USER=$POSTGRES_USER"
    echo "POSTGRES_PASSWORD=$POSTGRES_PASSWORD"
    echo "POSTGRES_DB=$POSTGRES_DB"

    # Authentication
    read -p "Enable Authentication? (y/n): " enable_auth
    if [[ "$enable_auth" == "y" ]]; then
      echo "ENABLE_AUTH=1"
      AUTH_METHOD_DEFAULT=${DEFAULT_AUTH_METHOD:-dev}
      read -p "Authentication method (dev/auth0/shibboleth) [${AUTH_METHOD_DEFAULT}]: " auth_method_input
      DEFAULT_AUTH_METHOD=$(echo "${auth_method_input:-$AUTH_METHOD_DEFAULT}" | tr '[:upper:]' '[:lower:]' | xargs)
      if [[ "$DEFAULT_AUTH_METHOD" == "dummy" ]]; then
        DEFAULT_AUTH_METHOD="dev"
      fi
      if [[ "$DEFAULT_AUTH_METHOD" != "dev" && "$DEFAULT_AUTH_METHOD" != "auth0" && "$DEFAULT_AUTH_METHOD" != "shibboleth" ]]; then
        echo -e "${YELLOW}Unknown auth method '$DEFAULT_AUTH_METHOD'. Falling back to dev.${NC}"
        DEFAULT_AUTH_METHOD="dev"
      fi
      OUTPUT_AUTH_METHOD="$DEFAULT_AUTH_METHOD"
      if [[ "$OUTPUT_AUTH_METHOD" == "dev" ]]; then
        OUTPUT_AUTH_METHOD="dummy"
      fi
      echo "DEFAULT_AUTH_METHOD=$OUTPUT_AUTH_METHOD"

      if [[ "$DEFAULT_AUTH_METHOD" == "dev" ]]; then
        DEV_SECRET_VALUE="${FLASK_SECRET_KEY:-change-me-in-production}"
        echo "FLASK_SECRET_KEY=$DEV_SECRET_VALUE"
        echo "$(prompt_variable "LOGIN_REDIRECT_URL" false "/login_dev")"
      elif [[ "$DEFAULT_AUTH_METHOD" == "auth0" ]]; then
        echo "$(prompt_variable "FLASK_SECRET_KEY" true "change-me-in-production")"
        echo "$(prompt_variable "LOGIN_REDIRECT_URL" false "/login")"
        echo "$(prompt_variable "AUTH0_DOMAIN" false)"
        echo "$(prompt_variable "AUTH0_CLIENT_ID" false)"
        echo "$(prompt_variable "AUTH0_CLIENT_SECRET" true)"
        echo "$(prompt_variable "AUTH0_CALLBACK_URL" false)"
        echo "$(prompt_variable "AUTH0_AUDIENCE" false)"
        echo "$(prompt_variable "AUTH0_DB_CONNECTION" false)"
        echo "$(prompt_variable "AUTH0_PUBLIC_KEY_URI" false)"
      else
        echo "$(prompt_variable "FLASK_SECRET_KEY" true "change-me-in-production")"
        echo "$(prompt_variable "LOGIN_REDIRECT_URL" false "/login_sso")"
        echo "$(prompt_variable "SHIBBOLETH_LOGIN_URL" false)"
        echo "$(prompt_variable "SHIBBOLETH_LOGOUT_URL" false)"
      fi
    else
      echo "ENABLE_AUTH=0"
      echo "DEFAULT_AUTH_METHOD=dummy"
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
  local sqlite_override_file="docker-sqlite.override.yml"
  local deploy_override_file=".docker-deploy.${DEPLOYMENT_NAME}.override.yml"
  local use_local_mdv_image="0"
  local mdv_mount="mdv-data:/app/mdv"
  local escaped_mount=""

  if [ ! -f "$compose_file" ]; then
    echo "Downloading $compose_file..."
    if ! curl -O "$compose_url"; then
      echo -e "${RED} Error: Failed to download $compose_file."
      exit 1
    fi
  fi

  if ! sed -E -i.bak "s/\"[0-9]+:5055\"/\"${APP_PORT}:5055\"/" "$compose_file"; then
    echo -e "${RED} Error: Failed to set app port in $compose_file.${NC}"
    exit 1
  fi

  if [[ "$STORAGE_MODE" == "host" ]]; then
    mkdir -p "$HOST_DATA_PATH"
    mdv_mount="${HOST_DATA_PATH}:/app/mdv"
  fi
  escaped_mount=${mdv_mount//\\/\\\\}
  escaped_mount=${escaped_mount//&/\\&}
  escaped_mount=${escaped_mount//|/\\|}
  escaped_mount=${escaped_mount//\"/\\\"}
  if ! sed -E -i.bak "s|^[[:space:]]*-[[:space:]]*.+:/app/mdv[[:space:]]*$|      - \"${escaped_mount}\"|" "$compose_file"; then
    echo -e "${RED} Error: Failed to set data mount in $compose_file.${NC}"
    exit 1
  fi

  cat > "$deploy_override_file" <<EOF
services:
  mdv_app:
    env_file:
      - "${ENV_FILE}"
EOF

  if [[ "${POSTGRES_MODE:-dedicated}" == "reuse" && -n "${REUSE_DB_NETWORK:-}" ]]; then
    cat >> "$deploy_override_file" <<EOF
    networks:
      - reused_db_net
EOF
  fi

  if [[ "${DB_BACKEND:-sqlite}" == "postgres" && "${POSTGRES_MODE:-dedicated}" != "reuse" ]]; then
    cat >> "$deploy_override_file" <<EOF
  mdv_db:
    container_name: mdv-db
    env_file:
      - "${ENV_FILE}"
EOF
  fi

  if [[ "${POSTGRES_MODE:-dedicated}" == "reuse" && -n "${REUSE_DB_NETWORK:-}" ]]; then
    cat >> "$deploy_override_file" <<EOF
networks:
  reused_db_net:
    external: true
    name: ${REUSE_DB_NETWORK}
EOF
  fi

  if docker image inspect mdvadmin/mdv:stable >/dev/null 2>&1; then
    use_local_mdv_image="1"
    echo "Detected local image mdvadmin/mdv:stable. Skipping mdv_app pull."
  fi

  if [[ "$DEPLOY_MODE" == "replace" ]]; then
    echo "Stopping existing deployment '${DEPLOYMENT_NAME}'..."
    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" down --remove-orphans || true
    if [[ "$RESET_DATA_VOLUME" == "1" ]]; then
      docker volume rm "${DEPLOYMENT_NAME}_mdv-data" >/dev/null 2>&1 || true
    fi
    if [[ "$REMOVE_OLD_VOLUME" == "1" ]]; then
      docker volume rm "${DEPLOYMENT_NAME}_mdv-data" >/dev/null 2>&1 || true
    fi
    if [[ "$REMOVE_OLD_HOST_PATH" == "1" && -n "$PREVIOUS_HOST_DATA_PATH" ]]; then
      rm -rf "$PREVIOUS_HOST_DATA_PATH" >/dev/null 2>&1 || true
    fi
  fi

  if [[ "${DB_BACKEND:-sqlite}" == "sqlite" ]]; then
    if [ ! -f "$sqlite_override_file" ]; then
      echo -e "${RED} Error: Missing $sqlite_override_file in project root.${NC}"
      exit 1
    fi
    if [[ "$use_local_mdv_image" != "1" ]]; then
      echo "Pulling mdv_app image..."
      docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" -f "$sqlite_override_file" pull mdv_app
    fi

    echo "Starting Docker Compose (sqlite backend)..."
    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" -f "$sqlite_override_file" up -d --no-deps mdv_app || {
      echo -e "${RED} Error: Failed to start Docker Compose for sqlite backend."
      exit 1
    }

    # Ensure legacy postgres service is not left running from previous deploys.
    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" stop mdv_db >/dev/null 2>&1 || true
    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" rm -f mdv_db >/dev/null 2>&1 || true
    docker volume rm "${DEPLOYMENT_NAME}_postgres-data" >/dev/null 2>&1 || true
    if [[ "${PREV_DB_BACKEND:-}" == "postgres" && "${PREV_DB_HOST:-}" == "mdv-db" && -n "${PREV_DB_SCHEMA:-}" ]]; then
      PREV_DB_USER_SAFE="${PREV_DB_USER:-testuser}"
      PREV_DB_NAME_SAFE="${PREV_DB_NAME:-mdv}"
      echo "Cleaning previous postgres schema '${PREV_DB_SCHEMA}' from mdv-db..."
      docker exec mdv-db psql -U "$PREV_DB_USER_SAFE" -d "$PREV_DB_NAME_SAFE" -c "DROP SCHEMA IF EXISTS \"${PREV_DB_SCHEMA}\" CASCADE;" >/dev/null 2>&1 || true
      NON_SYSTEM_SCHEMA_COUNT=$(docker exec mdv-db psql -U "$PREV_DB_USER_SAFE" -d "$PREV_DB_NAME_SAFE" -t -A -c "SELECT COUNT(*) FROM information_schema.schemata WHERE schema_name NOT IN ('pg_catalog','information_schema','public') AND schema_name NOT LIKE 'pg_toast%' AND schema_name NOT LIKE 'pg_temp_%';" 2>/dev/null | tr -d '[:space:]')
      if [[ "${NON_SYSTEM_SCHEMA_COUNT:-1}" == "0" ]]; then
        echo "No non-system schemas remain in mdv-db. Stopping and removing mdv-db container."
        docker stop mdv-db >/dev/null 2>&1 || true
        docker rm mdv-db >/dev/null 2>&1 || true
      fi
    fi
  elif [[ "${POSTGRES_MODE:-dedicated}" == "reuse" ]]; then
    if [[ "$use_local_mdv_image" != "1" ]]; then
      echo "Pulling mdv_app image..."
      docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" pull mdv_app
    fi

    echo "Starting Docker Compose (postgres reuse mode, app only)..."
    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" up -d --no-deps mdv_app || {
      echo -e "${RED} Error: Failed to start MDV app in postgres reuse mode.${NC}"
      exit 1
    }

    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" stop mdv_db >/dev/null 2>&1 || true
    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" rm -f mdv_db >/dev/null 2>&1 || true
  else
    if [[ "$use_local_mdv_image" == "1" ]]; then
      echo "Pulling postgres image only (mdv_app is local)..."
      docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" pull mdv_db
    else
      echo "Pulling stable Docker images..."
      docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" pull
    fi

    echo "Starting Docker Compose..."
    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" up -d || {
      echo -e "${RED} Error: Failed to start Docker Compose."
      exit 1
    }
  fi

  echo -e "${GREEN}✓ Docker Compose started successfully!"
}

# Function to open the application in the browser
open_browser() {
  URL="http://localhost:${APP_PORT}"
  echo -e "${GREEN}✓ Application deployed successfully!"
  echo -e "${GREEN} Open $URL in your browser."
}

# Main Execution
echo "Starting MDV Deployment..."

check_docker
#check_existing_service
prompt_deployment_settings

# Create app_logs directory if it doesn't exist
#mkdir -p app_logs && echo "Ensured app_logs directory exists."

create_or_validate_env_file

DOCKER_COMPOSE_URL="https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/main/docker-compose.yml"
run_docker_compose "$DOCKER_COMPOSE_URL"

open_browser
