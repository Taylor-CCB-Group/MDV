#!/bin/bash

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
#!/bin/bash

# Function to check and update environment variables
create_or_validate_env_file() {
  env_file=".env"

  if [ ! -f "$env_file" ]; then
    zenity --info --text="$env_file does not exist. Creating it now..."
    touch "$env_file"
  fi

  # Source the existing .env file to retain values
  source "$env_file"

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


  # Write the new .env file
    {
    echo "# Flask Configuration"
    echo "FLASK_ENV=production"
    echo "PYTHONUNBUFFERED=1"

    echo "# Database Configuration"
    DB_USER=$(prompt_variable "DB_USER" false)
    DB_PASSWORD=$(prompt_variable "DB_PASSWORD" true)
    DB_NAME=$(prompt_variable "DB_NAME" false)
    DB_HOST=$(prompt_variable "DB_HOST" false "mdv_db")  # Default value for DB_HOST

    # Output Database Variables
    echo "DB_USER=$DB_USER"
    echo "DB_PASSWORD=$DB_PASSWORD"
    echo "DB_NAME=$DB_NAME"
    echo "DB_HOST=$DB_HOST"

    # PostgreSQL should match DB values
    echo "POSTGRES_USER=$DB_USER"
    echo "POSTGRES_PASSWORD=$DB_PASSWORD"
    echo "POSTGRES_DB=$DB_NAME"

    # Authentication
    if zenity --question --title="Enable Authentication?" --text="Enable authentication?"; then
        echo "ENABLE_AUTH=1"
        FLASK_SECRET_KEY=$(prompt_variable "FLASK_SECRET_KEY" true)
        LOGIN_REDIRECT_URL=$(prompt_variable "LOGIN_REDIRECT_URL" false)
        AUTH0_DOMAIN=$(prompt_variable "AUTH0_DOMAIN" false)
        AUTH0_CLIENT_ID=$(prompt_variable "AUTH0_CLIENT_ID" false)
        AUTH0_CLIENT_SECRET=$(prompt_variable "AUTH0_CLIENT_SECRET" true)
        AUTH0_CALLBACK_URL=$(prompt_variable "AUTH0_CALLBACK_URL" false)
        AUTH0_AUDIENCE=$(prompt_variable "AUTH0_AUDIENCE" false)
        AUTH0_PUBLIC_KEY_URI=$(prompt_variable "AUTH0_PUBLIC_KEY_URI" false)

        # Shibboleth (only if authentication is enabled)
        SHIBBOLETH_LOGIN_URL=$(prompt_variable "SHIBBOLETH_LOGIN_URL" false)
        SHIBBOLETH_LOGOUT_URL=$(prompt_variable "SHIBBOLETH_LOGOUT_URL" false)

        # Output Authentication Variables
        echo "FLASK_SECRET_KEY=$FLASK_SECRET_KEY"
        echo "LOGIN_REDIRECT_URL=$LOGIN_REDIRECT_URL"
        echo "AUTH0_DOMAIN=$AUTH0_DOMAIN"
        echo "AUTH0_CLIENT_ID=$AUTH0_CLIENT_ID"
        echo "AUTH0_CLIENT_SECRET=$AUTH0_CLIENT_SECRET"
        echo "AUTH0_CALLBACK_URL=$AUTH0_CALLBACK_URL"
        echo "AUTH0_AUDIENCE=$AUTH0_AUDIENCE"
        echo "AUTH0_PUBLIC_KEY_URI=$AUTH0_PUBLIC_KEY_URI"
        echo "SHIBBOLETH_LOGIN_URL=$SHIBBOLETH_LOGIN_URL"
        echo "SHIBBOLETH_LOGOUT_URL=$SHIBBOLETH_LOGOUT_URL"
    else
        echo "ENABLE_AUTH=0"
    fi
    } > "$env_file"

  zenity --info --text="Environment variables updated and saved to $env_file."
}



# Function to download and run docker-compose
run_docker_compose() {
  local compose_url=$1
  local compose_file=$(basename "$compose_url")

  if [ ! -f "$compose_file" ]; then
    zenity --info --text="Downloading $compose_file..."
    if ! curl -O "$compose_url"; then
      zenity --error --text="Error: Failed to download $compose_file."
      exit 1
    fi
  fi

  zenity --info --text="Starting Docker Compose..."
  if ! docker compose -f "$compose_file" up -d; then
    zenity --error --text="Error: Failed to start Docker Compose."
    exit 1
  fi

  zenity --info --text="Docker Compose started successfully!"
}

# Function to open the browser
open_browser() {
  URL="https://localhost:5055"

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

# Create app_logs directory if it doesn't exist
if [ ! -d "app_logs" ]; then
    mkdir app_logs
    zenity --info --title="Directory Created" --text="The app_logs directory has been created."
else
    zenity --info --title="Directory Exists" --text="The app_logs directory already exists."
fi

create_or_validate_env_file

DOCKER_COMPOSE_URL="https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/auth_rbac/docker-secrets-local.yml"
run_docker_compose "$DOCKER_COMPOSE_URL"

zenity --info --title="Deployment Complete" --text="MDV application deployed successfully!\nClick OK to proceed."

open_browser

