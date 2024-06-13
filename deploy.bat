@echo off

REM Function to check if Docker is installed
:check_docker_installed
docker -v >nul 2>&1
if errorlevel 1 (
    echo Error: Docker is not installed. Please install Docker before running this script.
    exit /b 1
)

REM Function to check if Docker daemon is running
:check_docker_daemon
docker info >nul 2>&1
if errorlevel 1 (
    echo Error: Docker daemon is not running. Please start Docker before running this script.
    exit /b 1
)


REM Function to download and run docker-compose
:run_docker_compose
set PROD_COMPOSE_URL=https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/jh-dev/docker-compose.yml
echo Downloading docker-compose.yml from %PROD_COMPOSE_URL%...
curl -fsSL -o docker-compose.yml %PROD_COMPOSE_URL%
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

echo Welcome to the MDV application deployment script!

REM Call function to check if Docker is installed
call :check_docker_installed

REM Call function to check Docker daemon
call :check_docker_daemon

REM Call function to download and run docker-compose
call :run_docker_compose

echo MDV application deployment completed successfully!
exit /b 0
