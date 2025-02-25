#!/bin/bash

# Function - Check if Docker is installed
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

  if [ -f "$env_file" ]; then
      echo "$env_file already exists. Checking for required variables..."
      source "$env_file"

      # Check for database variables and prompt if missing
      if [ -z "$DB_USER" ]; then
          read -p "Please enter your database username: " DB_USER
      fi
      if [ -z "$DB_PASSWORD" ]; then
          read -p "Please enter your database password: " DB_PASSWORD
      fi
      if [ -z "$DB_NAME" ]; then
          read -p "Please enter your database name: " DB_NAME
      fi
      if [ -z "$DB_HOST" ]; then
          DB_HOST="mdv_db"  # Default value
      fi

      # Check for PostgreSQL variables and prompt if missing
      if [ -z "$POSTGRES_USER" ]; then
          POSTGRES_USER="$DB_USER"  # Default to DB_USER
      fi
      if [ -z "$POSTGRES_PASSWORD" ]; then
          POSTGRES_PASSWORD="$DB_PASSWORD"  # Default to DB_PASSWORD
      fi
      if [ -z "$POSTGRES_DB" ]; then
          POSTGRES_DB="$DB_NAME"  # Default to DB_NAME
      fi

      # Check if ENABLE_AUTH is set
      if [ -z "$ENABLE_AUTH" ]; then
          read -p "Do you want to enable the authentication feature? (yes/no): " enable_auth
          if [[ "$enable_auth" =~ ^[Yy] ]]; then
              ENABLE_AUTH=1
              # Prompt for FLASK_SECRET_KEY
              read -p "Please enter your Flask secret key: " FLASK_SECRET_KEY
              
              # Prompt for Auth0 variables if they are not already set
              if [ -z "$AUTH0_DOMAIN" ]; then
                  read -p "Please enter your Auth0 domain: " AUTH0_DOMAIN
              fi
              if [ -z "$AUTH0_CLIENT_ID" ]; then
                  read -p "Please enter your Auth0 client ID: " AUTH0_CLIENT_ID
              fi
              if [ -z "$AUTH0_CLIENT_SECRET" ]; then
                  read -p "Please enter your Auth0 client secret: " AUTH0_CLIENT_SECRET
              fi
              if [ -z "$AUTH0_CALLBACK_URL" ]; then
                  read -p "Please enter your Auth0 callback URL: " AUTH0_CALLBACK_URL
              fi
              if [ -z "$AUTH0_AUDIENCE" ]; then
                  read -p "Please enter your Auth0 audience: " AUTH0_AUDIENCE
              fi
              if [ -z "$AUTH0_PUBLIC_KEY_URI" ]; then
                  read -p "Please enter your Auth0 public key URI: " AUTH0_PUBLIC_KEY_URI
              fi
          else
              ENABLE_AUTH=0
              echo "Auth0 will not be enabled. Skipping Auth0 variable prompts."
          fi
      else
          echo "Authentication feature already configured with ENABLE_AUTH=$ENABLE_AUTH."
      fi

      # Check for Shibboleth variables if authentication is enabled
      if [ "$ENABLE_AUTH" -eq 1 ]; then
          if [ -z "$SHIBBOLETH_LOGIN_URL" ]; then
              read -p "Please enter your Shibboleth login URL: " SHIBBOLETH_LOGIN_URL
          fi
          if [ -z "$SHIBBOLETH_LOGOUT_URL" ]; then
              read -p "Please enter your Shibboleth logout URL: " SHIBBOLETH_LOGOUT_URL
          fi
      fi

  else
      echo "$env_file does not exist. Creating it now..."
      touch "$env_file"

      # Prompt for database credentials
      read -p "Please enter your database username: " DB_USER
      read -p "Please enter your database password: " DB_PASSWORD
      read -p "Please enter your database name: " DB_NAME
      DB_HOST="mdv_db"  # Default value

      # Prompt for PostgreSQL credentials
      POSTGRES_USER="$DB_USER"  # Default to DB_USER
      POSTGRES_PASSWORD="$DB_PASSWORD"  # Default to DB_PASSWORD
      POSTGRES_DB="$DB_NAME"  # Default to DB_NAME

      # Prompt for authentication feature
      read -p "Do you want to enable the authentication feature? (yes/no): " enable_auth
      if [[ "$enable_auth" =~ ^[Yy] ]]; then
          ENABLE_AUTH=1
          read -p "Please enter your Flask secret key: " FLASK_SECRET_KEY
          read -p "Please enter your Auth0 domain: " AUTH0_DOMAIN
          read -p "Please enter your Auth0 client ID: " AUTH0_CLIENT_ID
          read -p "Please enter your Auth0 client secret: " AUTH0_CLIENT_SECRET
          read -p "Please enter your Auth0 callback URL: " AUTH0_CALLBACK_URL
          read -p "Please enter your Auth0 audience: " AUTH0_AUDIENCE
          read -p "Please enter your Auth0 public key URI: " AUTH0_PUBLIC_KEY_URI
      else
          ENABLE_AUTH=0
          echo "Auth0 will not be enabled. Skipping Auth0 variable prompts."
      fi

      # Prompt for Shibboleth variables if authentication is enabled
      if [ "$ENABLE_AUTH" -eq 1 ]; then
          read -p "Please enter your Shibboleth login URL: " SHIBBOLETH_LOGIN_URL
          read -p "Please enter your Shibboleth logout URL: " SHIBBOLETH_LOGOUT_URL
      fi
  fi

  # Write the environment variables to the .env file
  {
      echo "# Flask Configuration"
      echo "FLASK_ENV=production"
      echo "PYTHONUNBUFFERED=1"
      echo "FLASK_SECRET_KEY=$FLASK_SECRET_KEY"  # Include FLASK_SECRET_KEY only if provided

      echo "# Database Configuration"
      echo "DB_USER=$DB_USER"
      echo "DB_PASSWORD=$DB_PASSWORD"
      echo "DB_NAME=$DB_NAME"
      echo "DB_HOST=$DB_HOST"  # Added host

      echo "POSTGRES_USER=$POSTGRES_USER"
      echo "POSTGRES_PASSWORD=$POSTGRES_PASSWORD"
      echo "POSTGRES_DB=$POSTGRES_DB"

      echo "# Authentication & Security"
      echo "LOGIN_REDIRECT_URL=https://localhost:5055/login_dev"
      echo "ENABLE_AUTH=$ENABLE_AUTH"
      if [ "$ENABLE_AUTH" -eq 1 ]; then
          echo "AUTH0_DOMAIN=$AUTH0_DOMAIN"
          echo "AUTH0_CLIENT_ID=$AUTH0_CLIENT_ID"
          echo "AUTH0_CLIENT_SECRET=$AUTH0_CLIENT_SECRET"
          echo "AUTH0_CALLBACK_URL=$AUTH0_CALLBACK_URL"
          echo "AUTH0_AUDIENCE=$AUTH0_AUDIENCE"
          echo "AUTH0_PUBLIC_KEY_URI=$AUTH0_PUBLIC_KEY_URI"
      fi

      echo "# Shibboleth Configuration"
      echo "SHIBBOLETH_LOGIN_URL=$SHIBBOLETH_LOGIN_URL"
      echo "SHIBBOLETH_LOGOUT_URL=$SHIBBOLETH_LOGOUT_URL"
  } > "$env_file" || {
    echo "Error: Failed to create $env_file."
    exit 1
  }

  echo "Environment variables saved to $env_file."

}

# Function to download and run docker-compose
run_docker_compose() {
  local compose_url=$1
  local compose_file=$(basename $compose_url)

  # Check if the compose file already exists
  if [ -f "$compose_file" ]; then
    echo "$compose_file already exists. Skipping download."
  else
    echo "Downloading $compose_file from $compose_url..."
    if ! curl -O "$compose_url"; then
      echo "Error: Failed to download $compose_file from $compose_url."
      exit 1
    fi
  fi

  echo "Running docker-compose with $compose_file..."
  if ! docker compose -f $compose_file up -d; then
    echo "Error: Failed to run docker-compose."
    exit 1
  fi
}

# Main execution starts here
echo "Welcome to the MDV application deployment script!"

# Check if Docker is installed
check_docker_installed

# Check if Docker daemon is running
check_docker_daemon

# Create or validate the .env file
create_or_validate_env_file

# URL of your production docker-compose file
DOCKER_COMPOSE_URL="https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/auth_rbac/docker-secrets-local.yml"

# Download and run docker-compose for production
run_docker_compose $DOCKER_COMPOSE_URL

echo "MDV application deployment completed successfully!"

echo
echo "******  Open your web browser and go to https://localhost:5055 to access the MDV application  ******"
echo
