@echo off

echo Welcome to the MDV application deployment script!

REM Function to check if Docker is installed
:check_docker_installed
echo Checking if Docker is installed...
docker -v >nul 2>&1
if errorlevel 1 (
    echo Error: Docker is not installed. Please install Docker before running this script.
    exit /b 1
)

REM Function to check if Docker daemon is running
:check_docker_daemon
echo Checking if Docker daemon is running...
docker info >nul 2>&1
if errorlevel 1 (
    echo Error: Docker daemon is not running. Please start Docker before running this script.
    exit /b 1
)

REM Function to download and run docker-compose
:run_docker_compose
echo Setting DOCKER_COMPOSE_URL environment variable...
set DOCKER_COMPOSE_URL=https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/jh-dev/docker-local.yml
echo Downloading docker-compose.yml from %DOCKER_COMPOSE_URL%...
curl -fsSL -o docker-compose.yml %DOCKER_COMPOSE_URL%
if errorlevel 1 (
    echo Error: Failed to download docker-compose.yml. Check the URL and try again.
    exit /b 1
)

echo Running docker-compose with docker-compose.yml...
docker-compose -f docker-compose.yml up -d
if errorlevel 1 (
    echo Error: Failed to run docker-compose. Check the configuration and try again.
    exit /b 1
)

echo MDV application deployment completed successfully!

echo.
echo ******  Open your web browser and go to https://localhost:5055 to access the MDV application  ******
echo.

exit /b 0
