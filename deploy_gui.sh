#!/bin/bash
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

# Function to check if Zenity is installed
check_zenity_installed() {
  if ! command -v zenity &> /dev/null; then
    echo "Error: Zenity is not installed. Please install it using your package manager."
    exit 1
  fi
}

# Function to check if Docker is installed
check_docker_installed() {
  if ! command -v docker &> /dev/null; then
    zenity --error --text="Error: Docker is not installed. Please install Docker before running this script."
    exit 1
  fi
}

# Function to check if Docker daemon is running
check_docker_daemon() {
  if ! docker info &> /dev/null; then
    zenity --error --text="Error: Docker daemon is not running. Please start Docker before running this script."
    exit 1
  fi
}

# Function to create or validate the .env file using Zenity
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
      if [ -n "$mount_type" ]; then
        echo "${mount_type}|${mount_source}"
      fi
    fi
  }

  DEPLOYMENT_NAME=$(zenity --entry --title="Deployment Name" --text="Deployment name:" --entry-text="mdv")
  if [ -z "$DEPLOYMENT_NAME" ]; then
    DEPLOYMENT_NAME="mdv"
  fi

  if docker ps -aq --filter "label=com.docker.compose.project=${DEPLOYMENT_NAME}" | grep -q .; then
    DEPLOY_ACTION=$(zenity --list --radiolist \
      --title="Existing Deployment Found" \
      --text="Deployment '${DEPLOYMENT_NAME}' already exists. What would you like to do?" \
      --column="Select" --column="Action" \
      TRUE "redeploy" FALSE "rename" FALSE "exit")
    if [ -z "$DEPLOY_ACTION" ] || [ "$DEPLOY_ACTION" = "exit" ]; then
      exit 0
    elif [ "$DEPLOY_ACTION" = "rename" ]; then
      DEPLOYMENT_NAME=$(zenity --entry --title="New Deployment Name" --text="Enter a new deployment name:" --entry-text="${DEPLOYMENT_NAME}_new")
      if [ -z "$DEPLOYMENT_NAME" ]; then
        zenity --error --text="Deployment name cannot be empty."
        exit 1
      fi
      if docker ps -aq --filter "label=com.docker.compose.project=${DEPLOYMENT_NAME}" | grep -q .; then
        zenity --error --text="Deployment '${DEPLOYMENT_NAME}' already exists. Choose a different name."
        exit 1
      fi
      DEPLOY_MODE="new"
    else
      DEPLOY_MODE="replace"
    fi
  else
    DEPLOY_MODE="new"
    zenity --info --text="No existing deployment found for '${DEPLOYMENT_NAME}'. Using new mode."
  fi

  DEFAULT_STORAGE_MODE="create"
  PREVIOUS_STORAGE_MOUNT=$(get_existing_storage_mount)
  if [ -n "$PREVIOUS_STORAGE_MOUNT" ]; then
    PREVIOUS_STORAGE_TYPE=${PREVIOUS_STORAGE_MOUNT%%|*}
    PREVIOUS_STORAGE_SOURCE=${PREVIOUS_STORAGE_MOUNT#*|}
    if [ "$PREVIOUS_STORAGE_TYPE" = "bind" ]; then
      PREVIOUS_STORAGE_MODE="host"
      PREVIOUS_HOST_DATA_PATH="$PREVIOUS_STORAGE_SOURCE"
      DEFAULT_STORAGE_MODE="host"
    elif [ "$PREVIOUS_STORAGE_TYPE" = "volume" ]; then
      PREVIOUS_STORAGE_MODE="reuse"
      DEFAULT_STORAGE_MODE="reuse"
    fi
  elif [ "$DEPLOY_MODE" = "replace" ] && docker volume inspect "${DEPLOYMENT_NAME}_mdv-data" >/dev/null 2>&1; then
    DEFAULT_STORAGE_MODE="reuse"
    PREVIOUS_STORAGE_MODE="reuse"
  fi

  APP_PORT_DEFAULT=$(get_next_available_port)
  if [ "$DEPLOY_MODE" = "replace" ]; then
    EXISTING_APP_PORT=$(get_existing_deployment_port)
    if [ -n "$EXISTING_APP_PORT" ]; then
      APP_PORT_DEFAULT="$EXISTING_APP_PORT"
      zenity --info --text="Detected existing app port $APP_PORT_DEFAULT for '$DEPLOYMENT_NAME'."
    fi
  fi
  APP_PORT=$(zenity --entry --title="Application Port" --text="Host port for MDV app:" --entry-text="$APP_PORT_DEFAULT")
  if [ -z "$APP_PORT" ]; then
    APP_PORT="$APP_PORT_DEFAULT"
  fi

  ENV_FILE=".env.${DEPLOYMENT_NAME}"
}

# Function to check and update environment variables
create_or_validate_env_file() {
  env_file="$ENV_FILE"

  if [ ! -f "$env_file" ]; then
    zenity --info --text="$env_file does not exist. Creating it now..."
    touch "$env_file"
  fi

  # Source the existing .env file to retain values
  source "$env_file"
  PREV_DB_BACKEND="${DB_BACKEND:-}"
  PREV_DB_HOST="${DB_HOST:-}"
  PREV_DB_SCHEMA="${DB_SCHEMA:-}"
  PREV_DB_USER="${DB_USER:-}"
  PREV_DB_NAME="${DB_NAME:-}"

  # Function to ask the user if they want to overwrite an existing variable
  prompt_variable() {
    local var_name=$1
    local hide_value=$2  # If "true", do not display the current value (for passwords, secrets)
    local default_value=$3  # Default value if variable is empty
    local current_value=${!var_name:-$default_value}  # Use existing value or default

    if [ -n "${!var_name}" ]; then
        # If the variable is already set, ask if the user wants to overwrite it
        if zenity --question --title="Overwrite $var_name?" --text="$var_name is already set.\nDo you want to overwrite it?"; then
            if [ "$hide_value" == "true" ]; then
                new_value=$(zenity --password --title="Enter new value for $var_name" --text="Enter a new value:")
            else
                new_value=$(zenity --entry --title="Enter new value for $var_name" --text="Current value: $current_value" --entry-text="$current_value")
            fi
        else
            new_value="${!var_name}"  # Keep the current value if not overwriting
        fi
    else
        # If there is no existing value, prompt for a new value
        if [ "$hide_value" == "true" ]; then
            new_value=$(zenity --password --title="Enter value for $var_name" --text="No existing value found. Enter a new value:")
        else
            new_value=$(zenity --entry --title="Enter value for $var_name" --text="No existing value found. Enter a new value:" --entry-text="$default_value")
        fi
    fi

    # Return only the new value
    echo "$new_value"
}

  ensure_shared_postgres() {
    if docker container inspect "$REUSE_DB_CONTAINER" >/dev/null 2>&1; then
      is_running=$(docker inspect -f '{{.State.Running}}' "$REUSE_DB_CONTAINER" 2>/dev/null || echo "false")
      if [ "$is_running" != "true" ]; then
        zenity --info --text="Shared postgres container '$REUSE_DB_CONTAINER' exists but is stopped. Starting it..."
        docker start "$REUSE_DB_CONTAINER" >/dev/null
      fi
      REUSE_DB_NETWORK=$(docker inspect "$REUSE_DB_CONTAINER" --format '{{range $k, $v := .NetworkSettings.Networks}}{{$k}}{{end}}')
      return 0
    fi

    if ! zenity --question --title="Create Shared Postgres" --text="Shared postgres container '$REUSE_DB_CONTAINER' not found.\nCreate it now?"; then
      zenity --error --text="Cannot continue in reuse mode without a shared postgres container."
      exit 1
    fi

    SHARED_DB_HOST_PATH=$(zenity --entry --title="Shared Postgres Data Folder" --text="Host folder for shared postgres data:" --entry-text="$SHARED_DB_HOST_PATH")
    if [ -z "$SHARED_DB_HOST_PATH" ]; then
      SHARED_DB_HOST_PATH="./mdv-data/postgres-shared"
    fi
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
      zenity --error --text="Failed to create shared postgres container '$REUSE_DB_CONTAINER'."
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
    zenity --error --text="Shared postgres container '$REUSE_DB_CONTAINER' did not become ready in time."
    exit 1
  }


  DB_BACKEND=$(zenity --list --radiolist \
    --title="Select Database Backend" \
    --text="Choose the database backend (default: sqlite)" \
    --column="Select" --column="Backend" \
    TRUE "sqlite" FALSE "postgres")
  if [ -z "$DB_BACKEND" ]; then
    DB_BACKEND="sqlite"
  fi

  SQLITE_DB_PATH=$(prompt_variable "SQLITE_DB_PATH" false "/app/mdv/mdv.sqlite3")

  if [ "$DB_BACKEND" = "postgres" ]; then
    DB_SCHEMA="$DEPLOYMENT_NAME"
    if docker container inspect "$REUSE_DB_CONTAINER" >/dev/null 2>&1; then
      REUSE_DB_NETWORK=$(docker inspect "$REUSE_DB_CONTAINER" --format '{{range $k, $v := .NetworkSettings.Networks}}{{$k}}{{end}}')
    fi
    EXISTING_POSTGRES_COUNT=0
    if docker container inspect "$REUSE_DB_CONTAINER" >/dev/null 2>&1; then
      EXISTING_POSTGRES_COUNT=1
    fi
    POSTGRES_MODE="reuse"
    zenity --info --text="Using shared postgres deployment mode (reuse)."
    DB_USER=$(prompt_variable "DB_USER" false "testuser")
    DB_PASSWORD=$(prompt_variable "DB_PASSWORD" true "testpass")
    DB_NAME=$(prompt_variable "DB_NAME" false "$DB_NAME")
    DB_HOST="$REUSE_DB_CONTAINER"
    ensure_shared_postgres
    DB_HOST=$(prompt_variable "DB_HOST" false "$DB_HOST")
  else
    DB_USER="${DB_USER:-testuser}"
    DB_PASSWORD="${DB_PASSWORD:-testpass}"
    DB_NAME="${DB_NAME:-mdv}"
    DB_HOST="mdv-db"
    DB_SCHEMA=""
  fi

  if [ "$DB_BACKEND" = "postgres" ]; then
    STORAGE_MODE="reuse"
    zenity --info --text="App data storage prompt skipped for postgres backend."
    if [ "${PREV_DB_BACKEND:-}" = "sqlite" ]; then
      if [ "$PREVIOUS_STORAGE_MODE" = "host" ] && [ -n "$PREVIOUS_HOST_DATA_PATH" ]; then
        if zenity --question --title="Remove Old Host Folder" --text="Switching sqlite->postgres.\nRemove old host folder '${PREVIOUS_HOST_DATA_PATH}' and its data?"; then
          REMOVE_OLD_HOST_PATH="1"
        fi
      fi
      if docker volume inspect "${DEPLOYMENT_NAME}_mdv-data" >/dev/null 2>&1; then
        if zenity --question --title="Remove Old App Volume" --text="Switching sqlite->postgres.\nRemove old app volume '${DEPLOYMENT_NAME}_mdv-data'?"; then
          REMOVE_OLD_VOLUME="1"
        fi
      fi
    fi
  else
    if [ "$DEPLOY_MODE" = "replace" ]; then
      if [ "$DEFAULT_STORAGE_MODE" = "host" ]; then
        zenity --info --text="Detected existing host-mapped data folder. Defaulting storage mode to host."
      elif [ "$DEFAULT_STORAGE_MODE" = "reuse" ]; then
        zenity --info --text="Detected existing volume-based data mount. Defaulting storage mode to reuse."
      fi
    fi
    if [ "$DEFAULT_STORAGE_MODE" = "reuse" ]; then
      STORAGE_MODE=$(zenity --list --radiolist --title="Data Storage Mode" --text="How should /app/mdv data be mounted?" --column="Select" --column="Mode" FALSE "create" TRUE "reuse" FALSE "host")
    else
      STORAGE_MODE=$(zenity --list --radiolist --title="Data Storage Mode" --text="How should /app/mdv data be mounted?" --column="Select" --column="Mode" TRUE "create" FALSE "reuse" FALSE "host")
    fi
    if [ -z "$STORAGE_MODE" ]; then
      STORAGE_MODE="$DEFAULT_STORAGE_MODE"
    fi
    if [ "$STORAGE_MODE" = "host" ]; then
      HOST_PATH_DEFAULT="./mdv-data/${DEPLOYMENT_NAME}"
      if [ "$PREVIOUS_STORAGE_MODE" = "host" ] && [ -n "$PREVIOUS_HOST_DATA_PATH" ]; then
        HOST_PATH_DEFAULT="$PREVIOUS_HOST_DATA_PATH"
      fi
      HOST_DATA_PATH=$(zenity --entry --title="Host Data Folder" --text="Host folder to mount at /app/mdv" --entry-text="$HOST_PATH_DEFAULT")
      if [ -z "$HOST_DATA_PATH" ]; then
        HOST_DATA_PATH="$HOST_PATH_DEFAULT"
      fi
      if [ "$PREVIOUS_STORAGE_MODE" = "reuse" ]; then
        if zenity --question --title="Remove Old Volume" --text="Switching from named volume to host folder.\nRemove old volume '${DEPLOYMENT_NAME}_mdv-data'?"; then
          REMOVE_OLD_VOLUME="1"
        fi
      fi
    elif [ "$STORAGE_MODE" = "create" ] && [ "$PREVIOUS_STORAGE_MODE" = "host" ] && [ -n "$PREVIOUS_HOST_DATA_PATH" ]; then
      if zenity --question --title="Remove Old Host Folder" --text="Switching from host folder to named volume.\nRemove old host folder '${PREVIOUS_HOST_DATA_PATH}' and its data?"; then
        REMOVE_OLD_HOST_PATH="1"
      fi
    fi
    if [ "$DEPLOY_MODE" = "replace" ] && [ "$STORAGE_MODE" = "create" ] && docker volume inspect "${DEPLOYMENT_NAME}_mdv-data" >/dev/null 2>&1; then
      if zenity --question --title="Existing Data Volume" --text="Existing data volume '${DEPLOYMENT_NAME}_mdv-data' found.\nRemove existing data before redeploy?"; then
        RESET_DATA_VOLUME="1"
      else
        STORAGE_MODE="reuse"
        zenity --info --text="Keeping existing data volume (reusing existing data)."
      fi
    fi
  fi

  if [ "$STORAGE_MODE" = "host" ]; then
    if [ "$PREVIOUS_STORAGE_MODE" = "host" ] && [ -n "$PREVIOUS_HOST_DATA_PATH" ] && [ "$PREVIOUS_HOST_DATA_PATH" != "$HOST_DATA_PATH" ]; then
      if zenity --question --title="Remove Old Host Folder" --text="Previous host folder differs:\n${PREVIOUS_HOST_DATA_PATH}\n\nRemove old folder and its data?"; then
        REMOVE_OLD_HOST_PATH="1"
      fi
    fi
    if [ -e "$HOST_DATA_PATH" ] && [ ! -d "$HOST_DATA_PATH" ]; then
      zenity --error --text="Error: '$HOST_DATA_PATH' exists but is not a folder."
      exit 1
    fi
    if [ -d "$HOST_DATA_PATH" ] && [ -n "$(ls -A "$HOST_DATA_PATH" 2>/dev/null)" ]; then
      if zenity --question --title="Host Folder Not Empty" --text="Host folder '$HOST_DATA_PATH' is not empty.\nRemove existing data before deployment?"; then
        rm -rf "$HOST_DATA_PATH"/* "$HOST_DATA_PATH"/.[!.]* "$HOST_DATA_PATH"/..?* 2>/dev/null || true
      else
        zenity --info --text="Keeping existing host folder data."
      fi
    fi
  fi

  # Write the new .env file
    {
    echo "# Flask Configuration"
    echo "FLASK_ENV=production"
    echo "PYTHONUNBUFFERED=1"

    echo "# Database Configuration"
    echo "DB_BACKEND=$DB_BACKEND"
    echo "POSTGRES_MODE=$POSTGRES_MODE"
    echo "SQLITE_DB_PATH=$SQLITE_DB_PATH"

    # Output Database Variables
    echo "DB_USER=$DB_USER"
    echo "DB_PASSWORD=$DB_PASSWORD"
    echo "DB_NAME=$DB_NAME"
    echo "DB_HOST=$DB_HOST"
    echo "DB_SCHEMA=$DB_SCHEMA"

    # PostgreSQL should match DB values
    echo "POSTGRES_USER=$DB_USER"
    echo "POSTGRES_PASSWORD=$DB_PASSWORD"
    echo "POSTGRES_DB=$DB_NAME"

    # Authentication
    if zenity --question --title="Enable Authentication?" --text="Enable authentication?"; then
        echo "ENABLE_AUTH=1"
        DEFAULT_AUTH_METHOD=$(zenity --list --radiolist \
          --title="Authentication Method" \
          --text="Select authentication provider" \
          --column="Select" --column="Method" \
          TRUE "dev" FALSE "auth0" FALSE "shibboleth")
        if [ -z "$DEFAULT_AUTH_METHOD" ]; then
            DEFAULT_AUTH_METHOD="dev"
        fi

        if [ "$DEFAULT_AUTH_METHOD" = "dummy" ]; then
            DEFAULT_AUTH_METHOD="dev"
        fi

        if [ "$DEFAULT_AUTH_METHOD" = "dev" ]; then
            FLASK_SECRET_KEY="${FLASK_SECRET_KEY:-change-me-in-production}"
            LOGIN_REDIRECT_URL=$(prompt_variable "LOGIN_REDIRECT_URL" false "/login_dev")
        elif [ "$DEFAULT_AUTH_METHOD" = "auth0" ]; then
            FLASK_SECRET_KEY=$(prompt_variable "FLASK_SECRET_KEY" true "change-me-in-production")
            LOGIN_REDIRECT_URL=$(prompt_variable "LOGIN_REDIRECT_URL" false "/login")
            AUTH0_DOMAIN=$(prompt_variable "AUTH0_DOMAIN" false)
            AUTH0_CLIENT_ID=$(prompt_variable "AUTH0_CLIENT_ID" false)
            AUTH0_CLIENT_SECRET=$(prompt_variable "AUTH0_CLIENT_SECRET" true)
            AUTH0_CALLBACK_URL=$(prompt_variable "AUTH0_CALLBACK_URL" false)
            AUTH0_AUDIENCE=$(prompt_variable "AUTH0_AUDIENCE" false)
            AUTH0_DB_CONNECTION=$(prompt_variable "AUTH0_DB_CONNECTION" false)
            AUTH0_PUBLIC_KEY_URI=$(prompt_variable "AUTH0_PUBLIC_KEY_URI" false)
        else
            FLASK_SECRET_KEY=$(prompt_variable "FLASK_SECRET_KEY" true "change-me-in-production")
            LOGIN_REDIRECT_URL=$(prompt_variable "LOGIN_REDIRECT_URL" false "/login_sso")
            SHIBBOLETH_LOGIN_URL=$(prompt_variable "SHIBBOLETH_LOGIN_URL" false)
            SHIBBOLETH_LOGOUT_URL=$(prompt_variable "SHIBBOLETH_LOGOUT_URL" false)
        fi

        # Output Authentication Variables (provider-specific keys are emitted only when configured)
        OUTPUT_AUTH_METHOD="$DEFAULT_AUTH_METHOD"
        if [ "$OUTPUT_AUTH_METHOD" = "dev" ]; then
            OUTPUT_AUTH_METHOD="dummy"
        fi
        echo "DEFAULT_AUTH_METHOD=$OUTPUT_AUTH_METHOD"
        echo "FLASK_SECRET_KEY=$FLASK_SECRET_KEY"
        echo "LOGIN_REDIRECT_URL=$LOGIN_REDIRECT_URL"
        if [ "$DEFAULT_AUTH_METHOD" = "auth0" ]; then
            echo "AUTH0_DOMAIN=$AUTH0_DOMAIN"
            echo "AUTH0_CLIENT_ID=$AUTH0_CLIENT_ID"
            echo "AUTH0_CLIENT_SECRET=$AUTH0_CLIENT_SECRET"
            echo "AUTH0_CALLBACK_URL=$AUTH0_CALLBACK_URL"
            echo "AUTH0_AUDIENCE=$AUTH0_AUDIENCE"
            echo "AUTH0_DB_CONNECTION=$AUTH0_DB_CONNECTION"
            echo "AUTH0_PUBLIC_KEY_URI=$AUTH0_PUBLIC_KEY_URI"
        elif [ "$DEFAULT_AUTH_METHOD" = "shibboleth" ]; then
            echo "SHIBBOLETH_LOGIN_URL=$SHIBBOLETH_LOGIN_URL"
            echo "SHIBBOLETH_LOGOUT_URL=$SHIBBOLETH_LOGOUT_URL"
        fi
    else
        echo "ENABLE_AUTH=0"
        echo "DEFAULT_AUTH_METHOD=dummy"
    fi
    } > "$env_file"

  zenity --info --text="Environment variables updated and saved to $env_file."
}



# Function to download and run docker-compose
run_docker_compose() {
  local compose_url=$1
  local compose_file=$(basename "$compose_url")
  local sqlite_override_file="docker-sqlite.override.yml"
  local deploy_override_file=".docker-deploy.${DEPLOYMENT_NAME}.override.yml"
  local use_local_mdv_image="0"
  local mdv_mount="mdv-data:/app/mdv"
  local escaped_mount=""

  if [ ! -f "$compose_file" ]; then
    zenity --info --text="Downloading $compose_file..."
    if ! curl -O "$compose_url"; then
      zenity --error --text="Error: Failed to download $compose_file."
      exit 1
    fi
  fi

  if ! sed -E -i.bak "s/\"[0-9]+:5055\"/\"${APP_PORT}:5055\"/" "$compose_file"; then
    zenity --error --text="Error: Failed to set app port in $compose_file."
    exit 1
  fi

  if [ "$STORAGE_MODE" = "host" ]; then
    mkdir -p "$HOST_DATA_PATH"
    mdv_mount="${HOST_DATA_PATH}:/app/mdv"
  fi
  escaped_mount=${mdv_mount//\\/\\\\}
  escaped_mount=${escaped_mount//&/\\&}
  escaped_mount=${escaped_mount//|/\\|}
  escaped_mount=${escaped_mount//\"/\\\"}
  if ! sed -E -i.bak "s|^[[:space:]]*-[[:space:]]*.+:/app/mdv[[:space:]]*$|      - \"${escaped_mount}\"|" "$compose_file"; then
    zenity --error --text="Error: Failed to set data mount in $compose_file."
    exit 1
  fi

  cat > "$deploy_override_file" <<EOF
services:
  mdv_app:
    env_file:
      - "${ENV_FILE}"
EOF

  if [ "${POSTGRES_MODE:-dedicated}" = "reuse" ] && [ -n "${REUSE_DB_NETWORK:-}" ]; then
    cat >> "$deploy_override_file" <<EOF
    networks:
      - reused_db_net
EOF
  fi

  if [ "${DB_BACKEND:-sqlite}" = "postgres" ] && [ "${POSTGRES_MODE:-dedicated}" != "reuse" ]; then
    cat >> "$deploy_override_file" <<EOF
  mdv_db:
    container_name: mdv-db
    env_file:
      - "${ENV_FILE}"
EOF
  fi

  if [ "${POSTGRES_MODE:-dedicated}" = "reuse" ] && [ -n "${REUSE_DB_NETWORK:-}" ]; then
    cat >> "$deploy_override_file" <<EOF
networks:
  reused_db_net:
    external: true
    name: ${REUSE_DB_NETWORK}
EOF
  fi

  if docker image inspect mdvadmin/mdv:stable >/dev/null 2>&1; then
    use_local_mdv_image="1"
    zenity --info --text="Detected local image mdvadmin/mdv:stable. Skipping mdv_app pull."
  fi

  if [ "$DEPLOY_MODE" = "replace" ]; then
    zenity --info --text="Stopping existing deployment '${DEPLOYMENT_NAME}'..."
    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" down --remove-orphans || true
    if [ "$RESET_DATA_VOLUME" = "1" ]; then
      docker volume rm "${DEPLOYMENT_NAME}_mdv-data" >/dev/null 2>&1 || true
    fi
    if [ "$REMOVE_OLD_VOLUME" = "1" ]; then
      docker volume rm "${DEPLOYMENT_NAME}_mdv-data" >/dev/null 2>&1 || true
    fi
    if [ "$REMOVE_OLD_HOST_PATH" = "1" ] && [ -n "$PREVIOUS_HOST_DATA_PATH" ]; then
      rm -rf "$PREVIOUS_HOST_DATA_PATH" >/dev/null 2>&1 || true
    fi
  fi

  if [ "${DB_BACKEND:-sqlite}" = "sqlite" ]; then
    if [ ! -f "$sqlite_override_file" ]; then
      zenity --error --text="Missing $sqlite_override_file. Please ensure it exists in the project root."
      exit 1
    fi

    if [ "$use_local_mdv_image" != "1" ]; then
      zenity --info --text="Pulling mdv_app image..."
      if ! docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" -f "$sqlite_override_file" pull mdv_app; then
        zenity --error --text="Error: Failed to pull mdv_app image."
        exit 1
      fi
    fi

    zenity --info --text="Starting Docker Compose (sqlite backend)..."
    if ! docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" -f "$sqlite_override_file" up -d --no-deps mdv_app; then
      zenity --error --text="Error: Failed to start Docker Compose for sqlite backend."
      exit 1
    fi

    # Ensure postgres service is not left running from previous deploys.
    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" stop mdv_db >/dev/null 2>&1 || true
    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" rm -f mdv_db >/dev/null 2>&1 || true
    docker volume rm "${DEPLOYMENT_NAME}_postgres-data" >/dev/null 2>&1 || true
    if [ "${PREV_DB_BACKEND:-}" = "postgres" ] && [ "${PREV_DB_HOST:-}" = "mdv-db" ] && [ -n "${PREV_DB_SCHEMA:-}" ]; then
      PREV_DB_USER_SAFE="${PREV_DB_USER:-testuser}"
      PREV_DB_NAME_SAFE="${PREV_DB_NAME:-mdv}"
      docker exec mdv-db psql -U "$PREV_DB_USER_SAFE" -d "$PREV_DB_NAME_SAFE" -c "DROP SCHEMA IF EXISTS \"${PREV_DB_SCHEMA}\" CASCADE;" >/dev/null 2>&1 || true
      NON_SYSTEM_SCHEMA_COUNT=$(docker exec mdv-db psql -U "$PREV_DB_USER_SAFE" -d "$PREV_DB_NAME_SAFE" -t -A -c "SELECT COUNT(*) FROM information_schema.schemata WHERE schema_name NOT IN ('pg_catalog','information_schema','public') AND schema_name NOT LIKE 'pg_toast%' AND schema_name NOT LIKE 'pg_temp_%';" 2>/dev/null | tr -d '[:space:]')
      if [ "${NON_SYSTEM_SCHEMA_COUNT:-1}" = "0" ]; then
        docker stop mdv-db >/dev/null 2>&1 || true
        docker rm mdv-db >/dev/null 2>&1 || true
      fi
    fi
  elif [ "${POSTGRES_MODE:-dedicated}" = "reuse" ]; then
    if [ "$use_local_mdv_image" != "1" ]; then
      zenity --info --text="Pulling mdv_app image..."
      if ! docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" pull mdv_app; then
        zenity --error --text="Error: Failed to pull mdv_app image."
        exit 1
      fi
    fi

    zenity --info --text="Starting Docker Compose (postgres reuse mode, app only)..."
    if ! docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" up -d --no-deps mdv_app; then
      zenity --error --text="Error: Failed to start Docker Compose for postgres reuse mode."
      exit 1
    fi

    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" stop mdv_db >/dev/null 2>&1 || true
    docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" rm -f mdv_db >/dev/null 2>&1 || true
  else
    if [ "$use_local_mdv_image" = "1" ]; then
      zenity --info --text="Pulling postgres image only (mdv_app is local)..."
      if ! docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" pull mdv_db; then
        zenity --error --text="Error: Failed to pull postgres image."
        exit 1
      fi
    else
      zenity --info --text="Pulling stable Docker images..."
      if ! docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" pull; then
        zenity --error --text="Error: Failed to pull Docker images."
        exit 1
      fi
    fi

    zenity --info --text="Starting Docker Compose..."
    if ! docker compose -p "$DEPLOYMENT_NAME" --env-file "$ENV_FILE" -f "$compose_file" -f "$deploy_override_file" up -d; then
      zenity --error --text="Error: Failed to start Docker Compose."
      exit 1
    fi
  fi

  zenity --info --text="Docker Compose started successfully!"
}

# Function to open the browser
open_browser() {
  URL="http://localhost:${APP_PORT}"

  zenity --question --title="Open Application" --text="Do you want to open MDV in your browser now?"
  if [ $? -eq 0 ]; then
    if command -v xdg-open &> /dev/null; then
      xdg-open "$URL"  # Linux
    elif command -v open &> /dev/null; then
      open "$URL"  # macOS
    else
      zenity --error --text="Could not detect a way to open the browser. Please open it manually: $URL"
    fi
  fi
}

# Main Execution
zenity --info --title="MDV Deployment" --text="Welcome to the MDV application deployment script!"

check_zenity_installed
check_docker_installed
check_docker_daemon
prompt_deployment_settings

# Create app_logs directory if it doesn't exist
if [ ! -d "app_logs" ]; then
    mkdir app_logs
    zenity --info --title="Directory Created" --text="The app_logs directory has been created."
else
    zenity --info --title="Directory Exists" --text="The app_logs directory already exists."
fi

create_or_validate_env_file

DOCKER_COMPOSE_URL="https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/main/docker-compose.yml"
run_docker_compose "$DOCKER_COMPOSE_URL"

zenity --info --title="Deployment Complete" --text="MDV application deployed successfully!\nClick OK to proceed."

open_browser
