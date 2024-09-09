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

# Function to download and run docker-compose
run_docker_compose() {
  local compose_url=$1
  local compose_file=$(basename $compose_url)
  
  echo "Downloading $compose_file from $compose_url..."
  curl -O $compose_url  # Download the file
  
  echo "Running docker-compose with $compose_file..."
  docker-compose -f $compose_file up -d  # Run docker-compose with the downloaded file
}

# Main execution starts here
echo "Welcome to the MDV application deployment script!"

# Check if Docker is installed
check_docker_installed

# Check if Docker daemon is running
check_docker_daemon

# URL of your production docker-compose file
DOCKER_COMPOSE_URL="https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/jh-dev/docker-local.yml"

# Download and run docker-compose for production
run_docker_compose $DOCKER_COMPOSE_URL

echo "MDV application deployment completed successfully!"

echo
echo "******  Open your web browser and go to https://localhost:5055 to access the MDV application  ******"
echo
