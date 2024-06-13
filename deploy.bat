@echo off
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

echo Welcome to the MDV (Market Data Validation) application deployment script!

REM Check if Docker is installed
docker -v >nul 2>&1
if errorlevel 1 (
    echo Error: Docker is not installed. Please install Docker before running this script.
    exit /b 1
)

REM Call function to check Docker daemon
call :check_docker_daemon

REM Download and run docker-compose for production
echo Downloading docker-compose.yml from %PROD_COMPOSE_URL%...
curl -O %PROD_COMPOSE_URL%

echo Running docker-compose with docker-compose.yml...
docker-compose -f docker-compose.yml up -d

echo MDV application deployment completed successfully!
exit /b 0
