@echo off
setlocal EnableDelayedExpansion

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

echo.
set /p DEPLOYMENT_NAME=Deployment name [mdv]:
if "%DEPLOYMENT_NAME%"=="" set DEPLOYMENT_NAME=mdv
set HAS_EXISTING_DEPLOYMENT=0
for /f %%i in ('docker ps -aq --filter "label=com.docker.compose.project=%DEPLOYMENT_NAME%"') do set HAS_EXISTING_DEPLOYMENT=1
if "%HAS_EXISTING_DEPLOYMENT%"=="1" (
    set /p DEPLOY_ACTION=Deployment "%DEPLOYMENT_NAME%" exists. Action ^(redeploy/rename/exit^) [redeploy]:
    if /I "!DEPLOY_ACTION!"=="" set DEPLOY_ACTION=redeploy
    if /I "!DEPLOY_ACTION!"=="exit" exit /b 0
    if /I "!DEPLOY_ACTION!"=="rename" (
        set /p DEPLOYMENT_NAME=New deployment name:
        if "!DEPLOYMENT_NAME!"=="" (
            echo Error: deployment name cannot be empty.
            exit /b 1
        )
        set HAS_EXISTING_DEPLOYMENT=0
        for /f %%i in ('docker ps -aq --filter "label=com.docker.compose.project=!DEPLOYMENT_NAME!"') do set HAS_EXISTING_DEPLOYMENT=1
        if "!HAS_EXISTING_DEPLOYMENT!"=="1" (
            echo Error: deployment "!DEPLOYMENT_NAME!" already exists. Please choose a different name.
            exit /b 1
        )
        set DEPLOY_MODE=new
    ) else (
        set DEPLOY_MODE=replace
    )
) else (
    set DEPLOY_MODE=new
    echo No existing deployment found for "%DEPLOYMENT_NAME%". Using new mode.
)
set ENV_FILE=.env.%DEPLOYMENT_NAME%

set PREVIOUS_STORAGE_MODE=
set PREVIOUS_HOST_DATA_PATH=
if /I "%DEPLOY_MODE%"=="replace" (
    for /f "tokens=1,* delims=|" %%A in ('powershell -NoProfile -Command "$m = docker inspect '%DEPLOYMENT_NAME%-mdv_app-1' 2>$null | ConvertFrom-Json; if($m){ $mount = $m[0].Mounts | Where-Object {$_.Destination -eq '/app/mdv'} | Select-Object -First 1; if($mount){ if($mount.Type -eq 'bind'){ Write-Output ('host|' + $mount.Source) } elseif($mount.Type -eq 'volume'){ Write-Output 'reuse|' } } }"') do (
        set PREVIOUS_STORAGE_MODE=%%A
        set PREVIOUS_HOST_DATA_PATH=%%B
    )
)
set DEFAULT_STORAGE_MODE=create
if /I "%DEPLOY_MODE%"=="replace" (
    if /I "!PREVIOUS_STORAGE_MODE!"=="host" (
        set DEFAULT_STORAGE_MODE=host
    ) else if /I "!PREVIOUS_STORAGE_MODE!"=="reuse" (
        set DEFAULT_STORAGE_MODE=reuse
    ) else (
        docker volume inspect %DEPLOYMENT_NAME%_mdv-data >nul 2>&1
        if not errorlevel 1 set DEFAULT_STORAGE_MODE=reuse
    )
)
set STORAGE_MODE=
set HOST_DATA_PATH=
set RESET_DATA_VOLUME=0
set REMOVE_OLD_HOST_PATH=0
set REMOVE_OLD_VOLUME=0

set APP_PORT_DEFAULT=5055
for /f %%p in ('powershell -NoProfile -Command "$ports = docker ps --format '{{.Ports}}'; $p=5055; while($true){ if(-not($ports -match ('0.0.0.0:'+$p+'->') -or $ports -match ('\[\:\:\]:'+$p+'->'))){ break }; $p++ }; Write-Output $p"') do set APP_PORT_DEFAULT=%%p
if /I "%DEPLOY_MODE%"=="replace" (
    set EXISTING_APP_PORT=
    for /f %%p in ('powershell -NoProfile -Command "$line = docker port '%DEPLOYMENT_NAME%-mdv_app-1' 5055/tcp 2>$null | Select-Object -First 1; if($line -match ':(\d+)$'){ Write-Output $matches[1] }"') do (
        set APP_PORT_DEFAULT=%%p
        set EXISTING_APP_PORT=%%p
    )
    if not "!EXISTING_APP_PORT!"=="" echo Detected existing app port !APP_PORT_DEFAULT! for "%DEPLOYMENT_NAME%".
)
set /p APP_PORT=Application host port [%APP_PORT_DEFAULT%]:
if "%APP_PORT%"=="" set APP_PORT=%APP_PORT_DEFAULT%

set PREV_DB_BACKEND=
set PREV_DB_HOST=
set PREV_DB_SCHEMA=
set PREV_DB_USER=
set PREV_DB_PASSWORD=
set PREV_DB_NAME=
if exist "%ENV_FILE%" (
    for /f "usebackq tokens=1,* delims==" %%A in ("%ENV_FILE%") do (
        if /I "%%A"=="DB_BACKEND" set PREV_DB_BACKEND=%%B
        if /I "%%A"=="DB_HOST" set PREV_DB_HOST=%%B
        if /I "%%A"=="DB_SCHEMA" set PREV_DB_SCHEMA=%%B
        if /I "%%A"=="DB_USER" set PREV_DB_USER=%%B
        if /I "%%A"=="DB_PASSWORD" set PREV_DB_PASSWORD=%%B
        if /I "%%A"=="DB_NAME" set PREV_DB_NAME=%%B
    )
)

REM Configure .env with database backend selection (default sqlite)
echo.
set /p DB_BACKEND=Choose database backend ^(sqlite/postgres^) [sqlite]:
if "%DB_BACKEND%"=="" set DB_BACKEND=sqlite
if /I not "%DB_BACKEND%"=="sqlite" if /I not "%DB_BACKEND%"=="postgres" (
    echo Invalid backend "%DB_BACKEND%". Falling back to sqlite.
    set DB_BACKEND=sqlite
)

set DB_HOST=mdv-db
if "%SQLITE_DB_PATH%"=="" set SQLITE_DB_PATH=/app/mdv/mdv.sqlite3
set POSTGRES_MODE=dedicated
set REUSE_DB_CONTAINER=mdv-db
set REUSE_DB_NETWORK=
set DB_SCHEMA=%DEPLOYMENT_NAME%
set SHARED_PG_ENV=.env.postgres-shared
set DB_NAME=mdv
set SHARED_DB_NETWORK=mdv-shared-db
set SHARED_DB_HOST_PATH=%CD%\mdv-data\postgres-shared

if /I "%DB_BACKEND%"=="postgres" (
    set EXISTING_POSTGRES=0
    docker container inspect mdv-db >nul 2>&1
    if not errorlevel 1 set EXISTING_POSTGRES=1
    if "!EXISTING_POSTGRES!"=="1" (
        for /f %%i in ('docker inspect mdv-db --format "{{range $k, $v := .NetworkSettings.Networks}}{{$k}}{{end}}"') do set REUSE_DB_NETWORK=%%i
    )
    set POSTGRES_MODE=reuse
    echo Using shared postgres deployment mode ^(reuse^).
    if "!EXISTING_POSTGRES!"=="1" (
        echo Detected running shared Postgres container ^(mdv-db^).
    )
    if /I "!POSTGRES_MODE!"=="reuse" (
        set DB_USER=
        set DB_PASSWORD=
        if exist "!SHARED_PG_ENV!" (
            for /f "usebackq tokens=1,* delims==" %%A in ("!SHARED_PG_ENV!") do (
                if /I "%%A"=="DB_USER" set DB_USER=%%B
                if /I "%%A"=="DB_PASSWORD" set DB_PASSWORD=%%B
            )
        )
        if "!DB_USER!"=="" if "!EXISTING_POSTGRES!"=="1" for /f %%i in ('powershell -NoProfile -Command "$envs=(docker inspect mdv-db --format '{{json .Config.Env}}' ^| ConvertFrom-Json); ($envs ^| ? {$_ -like 'POSTGRES_USER=*'} ^| %% {$_.Split('=')[1]})"') do set DB_USER=%%i
        if "!DB_PASSWORD!"=="" if "!EXISTING_POSTGRES!"=="1" for /f %%i in ('powershell -NoProfile -Command "$envs=(docker inspect mdv-db --format '{{json .Config.Env}}' ^| ConvertFrom-Json); ($envs ^| ? {$_ -like 'POSTGRES_PASSWORD=*'} ^| %% {$_.Split('=')[1]})"') do set DB_PASSWORD=%%i
        if "!DB_USER!"=="" set DB_USER=testuser
        if "!DB_PASSWORD!"=="" set DB_PASSWORD=testpass
        set /p DB_NAME=Enter DB_NAME [!DB_NAME!]:
        if "!DB_NAME!"=="" set DB_NAME=mdv

        if "!EXISTING_POSTGRES!"=="0" (
            set /p CREATE_SHARED_DB_INPUT=Shared postgres container "mdv-db" not found. Create it now? ^(y/n^) [y]:
            if /I "!CREATE_SHARED_DB_INPUT!"=="n" (
                echo Cannot continue in reuse mode without shared postgres.
                exit /b 1
            )
            set /p SHARED_DB_HOST_PATH=Shared Postgres host folder [!SHARED_DB_HOST_PATH!]:
            if "!SHARED_DB_HOST_PATH!"=="" set SHARED_DB_HOST_PATH=%CD%\mdv-data\postgres-shared
            docker network inspect !SHARED_DB_NETWORK! >nul 2>&1
            if errorlevel 1 docker network create !SHARED_DB_NETWORK! >nul 2>&1
            if not exist "!SHARED_DB_HOST_PATH!" mkdir "!SHARED_DB_HOST_PATH!"
            docker run -d --name mdv-db --network !SHARED_DB_NETWORK! -e POSTGRES_USER=!DB_USER! -e POSTGRES_PASSWORD=!DB_PASSWORD! -e POSTGRES_DB=!DB_NAME! -v "!SHARED_DB_HOST_PATH!:/var/lib/postgresql/data" postgres:16 >nul 2>&1
            if errorlevel 1 (
                echo Failed to create shared postgres container mdv-db.
                exit /b 1
            )
            set REUSE_DB_NETWORK=!SHARED_DB_NETWORK!
            set EXISTING_POSTGRES=1
        ) else (
            for /f %%i in ('docker inspect mdv-db --format "{{.State.Running}}"') do set SHARED_DB_RUNNING=%%i
            if /I not "!SHARED_DB_RUNNING!"=="true" (
                echo Shared postgres container mdv-db exists but is stopped. Starting it...
                docker start mdv-db >nul 2>&1
            )
        )
        if "!REUSE_DB_NETWORK!"=="" for /f %%i in ('docker inspect mdv-db --format "{{range $k, $v := .NetworkSettings.Networks}}{{$k}}{{end}}"') do set REUSE_DB_NETWORK=%%i
        set DEFAULT_DB_HOST=mdv-db
        set /p DB_HOST=Existing Postgres host [!DEFAULT_DB_HOST!]:
        if "!DB_HOST!"=="" set DB_HOST=!DEFAULT_DB_HOST!
        echo Reusing saved/existing Postgres credentials for convenience.
    )
) else (
    if "!DB_USER!"=="" set DB_USER=testuser
    if "!DB_PASSWORD!"=="" set DB_PASSWORD=testpass
    if "!DB_NAME!"=="" set DB_NAME=mdv
    set DB_SCHEMA=
)

if /I "%DB_BACKEND%"=="postgres" (
    set STORAGE_MODE=reuse
    echo App data storage is not prompted for postgres backend ^(shared postgres deployment^).
    if /I "!PREV_DB_BACKEND!"=="sqlite" (
        if /I "!PREVIOUS_STORAGE_MODE!"=="host" if /I not "!PREVIOUS_HOST_DATA_PATH!"=="" (
            set /p REMOVE_OLD_HOST_INPUT=Switching sqlite->postgres. Remove old host folder "!PREVIOUS_HOST_DATA_PATH!" and its data? ^(y/n^) [n]:
            if /I "!REMOVE_OLD_HOST_INPUT!"=="y" set REMOVE_OLD_HOST_PATH=1
        )
        docker volume inspect %DEPLOYMENT_NAME%_mdv-data >nul 2>&1
        if not errorlevel 1 (
            set /p REMOVE_OLD_VOLUME_INPUT=Switching sqlite->postgres. Remove old app volume "%DEPLOYMENT_NAME%_mdv-data"? ^(y/n^) [n]:
            if /I "!REMOVE_OLD_VOLUME_INPUT!"=="y" set REMOVE_OLD_VOLUME=1
        )
    )
) else (
    if /I "%DEPLOY_MODE%"=="replace" (
        if /I "!DEFAULT_STORAGE_MODE!"=="host" (
            echo Detected existing host-mapped data folder. Defaulting storage mode to host.
        ) else if /I "!DEFAULT_STORAGE_MODE!"=="reuse" (
            echo Detected existing volume-based data mount. Defaulting storage mode to reuse.
        )
    )
    set /p STORAGE_MODE=Data storage mode ^(create/reuse/host^) [%DEFAULT_STORAGE_MODE%]:
    if "!STORAGE_MODE!"=="" set STORAGE_MODE=!DEFAULT_STORAGE_MODE!
    if /I not "!STORAGE_MODE!"=="create" if /I not "!STORAGE_MODE!"=="reuse" if /I not "!STORAGE_MODE!"=="host" (
        echo Invalid storage mode "!STORAGE_MODE!". Falling back to !DEFAULT_STORAGE_MODE!.
        set STORAGE_MODE=!DEFAULT_STORAGE_MODE!
    )
    if /I "!STORAGE_MODE!"=="host" (
        set HOST_PATH_DEFAULT=.\mdv-data\%DEPLOYMENT_NAME%
        if /I "!PREVIOUS_STORAGE_MODE!"=="host" if not "!PREVIOUS_HOST_DATA_PATH!"=="" set HOST_PATH_DEFAULT=!PREVIOUS_HOST_DATA_PATH!
        set /p HOST_DATA_PATH=Host folder for /app/mdv [!HOST_PATH_DEFAULT!]:
        if "!HOST_DATA_PATH!"=="" set HOST_DATA_PATH=!HOST_PATH_DEFAULT!
        if /I "!PREVIOUS_STORAGE_MODE!"=="reuse" (
            set /p REMOVE_OLD_VOLUME_INPUT=Switching from named volume to host folder. Remove old volume "%DEPLOYMENT_NAME%_mdv-data"? ^(y/n^) [n]:
            if /I "!REMOVE_OLD_VOLUME_INPUT!"=="y" set REMOVE_OLD_VOLUME=1
        )
    ) else if /I "!STORAGE_MODE!"=="create" (
        if /I "!PREVIOUS_STORAGE_MODE!"=="host" if /I not "!PREVIOUS_HOST_DATA_PATH!"=="" (
            set /p REMOVE_OLD_HOST_INPUT=Switching from host folder to named volume. Remove old host folder "!PREVIOUS_HOST_DATA_PATH!" and its data? ^(y/n^) [n]:
            if /I "!REMOVE_OLD_HOST_INPUT!"=="y" set REMOVE_OLD_HOST_PATH=1
        )
    )
    if /I "%DEPLOY_MODE%"=="replace" if /I "!STORAGE_MODE!"=="create" (
        docker volume inspect %DEPLOYMENT_NAME%_mdv-data >nul 2>&1
        if not errorlevel 1 (
            set /p REMOVE_DATA_INPUT=Existing data volume "%DEPLOYMENT_NAME%_mdv-data" found. Remove existing data before redeploy? ^(y/n^) [n]:
            if /I "!REMOVE_DATA_INPUT!"=="y" (
                set RESET_DATA_VOLUME=1
            ) else (
                set STORAGE_MODE=reuse
                echo Keeping existing data volume ^(reusing existing data^).
            )
        )
    )
)

if /I "%STORAGE_MODE%"=="host" (
    if /I "!PREVIOUS_STORAGE_MODE!"=="host" if /I not "!PREVIOUS_HOST_DATA_PATH!"=="" if /I not "!PREVIOUS_HOST_DATA_PATH!"=="!HOST_DATA_PATH!" (
        set /p REMOVE_OLD_HOST_INPUT=Previous host folder "!PREVIOUS_HOST_DATA_PATH!" differs. Remove old folder and its data? ^(y/n^) [n]:
        if /I "!REMOVE_OLD_HOST_INPUT!"=="y" set REMOVE_OLD_HOST_PATH=1
    )
    if exist "!HOST_DATA_PATH!\." (
        set HOST_PATH_NOT_EMPTY=0
        for /f %%i in ('dir /b "!HOST_DATA_PATH!" 2^>nul') do set HOST_PATH_NOT_EMPTY=1
        if "!HOST_PATH_NOT_EMPTY!"=="1" (
            set /p CLEAR_HOST_INPUT=Host folder "!HOST_DATA_PATH!" is not empty. Remove existing data before deployment? ^(y/n^) [n]:
            if /I "!CLEAR_HOST_INPUT!"=="y" (
                del /f /q "!HOST_DATA_PATH!\*" >nul 2>&1
                for /d %%D in ("!HOST_DATA_PATH!\*") do rd /s /q "%%D"
            ) else (
                echo Keeping existing host folder data.
            )
        )
    ) else if exist "!HOST_DATA_PATH!" (
        echo Error: "!HOST_DATA_PATH!" exists but is not a folder.
        exit /b 1
    ) else (
        mkdir "!HOST_DATA_PATH!"
    )
)

set POSTGRES_USER=!DB_USER!
set POSTGRES_PASSWORD=!DB_PASSWORD!
set POSTGRES_DB=!DB_NAME!

if /I "!DB_BACKEND!"=="postgres" (
    (
        echo DB_USER=!DB_USER!
        echo DB_PASSWORD=!DB_PASSWORD!
    ) > "!SHARED_PG_ENV!"
)

echo.
set ENABLE_AUTH=0
set DEFAULT_AUTH_METHOD=dev
set DEFAULT_AUTH_METHOD_ENV=dummy
set FLASK_SECRET_KEY=
set LOGIN_REDIRECT_URL=
set AUTH0_DOMAIN=
set AUTH0_CLIENT_ID=
set AUTH0_CLIENT_SECRET=
set AUTH0_CALLBACK_URL=
set AUTH0_AUDIENCE=
set AUTH0_DB_CONNECTION=
set AUTH0_PUBLIC_KEY_URI=
set SHIBBOLETH_LOGIN_URL=
set SHIBBOLETH_LOGOUT_URL=

set /p ENABLE_AUTH_INPUT=Enable Authentication? ^(y/n^) [n]:
if /I "!ENABLE_AUTH_INPUT!"=="y" (
    set ENABLE_AUTH=1
    set /p DEFAULT_AUTH_METHOD=Authentication method ^(dev/auth0/shibboleth^) [dev]:
    if "!DEFAULT_AUTH_METHOD!"=="" set DEFAULT_AUTH_METHOD=dev
    if /I "!DEFAULT_AUTH_METHOD!"=="dummy" set DEFAULT_AUTH_METHOD=dev
    if /I not "!DEFAULT_AUTH_METHOD!"=="dev" if /I not "!DEFAULT_AUTH_METHOD!"=="auth0" if /I not "!DEFAULT_AUTH_METHOD!"=="shibboleth" (
        echo Invalid auth method "!DEFAULT_AUTH_METHOD!". Falling back to dev.
        set DEFAULT_AUTH_METHOD=dev
    )

    if /I "!DEFAULT_AUTH_METHOD!"=="dev" (
        set FLASK_SECRET_KEY=change-me-in-production
        set /p LOGIN_REDIRECT_URL=Enter LOGIN_REDIRECT_URL [/login_dev]:
        if "!LOGIN_REDIRECT_URL!"=="" set LOGIN_REDIRECT_URL=/login_dev
    ) else if /I "!DEFAULT_AUTH_METHOD!"=="auth0" (
        set /p FLASK_SECRET_KEY=Enter FLASK_SECRET_KEY [change-me-in-production]:
        if "!FLASK_SECRET_KEY!"=="" set FLASK_SECRET_KEY=change-me-in-production
        set /p LOGIN_REDIRECT_URL=Enter LOGIN_REDIRECT_URL [/login]:
        if "!LOGIN_REDIRECT_URL!"=="" set LOGIN_REDIRECT_URL=/login

        set /p AUTH0_DOMAIN=Enter AUTH0_DOMAIN:
        set /p AUTH0_CLIENT_ID=Enter AUTH0_CLIENT_ID:
        set /p AUTH0_CLIENT_SECRET=Enter AUTH0_CLIENT_SECRET:
        set /p AUTH0_CALLBACK_URL=Enter AUTH0_CALLBACK_URL:
        set /p AUTH0_AUDIENCE=Enter AUTH0_AUDIENCE:
        set /p AUTH0_DB_CONNECTION=Enter AUTH0_DB_CONNECTION:
        set /p AUTH0_PUBLIC_KEY_URI=Enter AUTH0_PUBLIC_KEY_URI:
    ) else (
        set /p FLASK_SECRET_KEY=Enter FLASK_SECRET_KEY [change-me-in-production]:
        if "!FLASK_SECRET_KEY!"=="" set FLASK_SECRET_KEY=change-me-in-production
        set /p LOGIN_REDIRECT_URL=Enter LOGIN_REDIRECT_URL [/login_sso]:
        if "!LOGIN_REDIRECT_URL!"=="" set LOGIN_REDIRECT_URL=/login_sso

        set /p SHIBBOLETH_LOGIN_URL=Enter SHIBBOLETH_LOGIN_URL:
        set /p SHIBBOLETH_LOGOUT_URL=Enter SHIBBOLETH_LOGOUT_URL:
    )
)

set DEFAULT_AUTH_METHOD_ENV=!DEFAULT_AUTH_METHOD!
if /I "!DEFAULT_AUTH_METHOD_ENV!"=="dev" set DEFAULT_AUTH_METHOD_ENV=dummy

echo Writing .env file...
(
    echo # Flask Configuration
    echo FLASK_ENV=production
    echo PYTHONUNBUFFERED=1
    echo.
    echo # Database Configuration
    echo DB_BACKEND=!DB_BACKEND!
    echo POSTGRES_MODE=!POSTGRES_MODE!
    echo SQLITE_DB_PATH=!SQLITE_DB_PATH!
    echo DB_USER=!DB_USER!
    echo DB_PASSWORD=!DB_PASSWORD!
    echo DB_NAME=!DB_NAME!
    echo DB_HOST=!DB_HOST!
    echo DB_SCHEMA=!DB_SCHEMA!
    echo.
    echo POSTGRES_USER=!POSTGRES_USER!
    echo POSTGRES_PASSWORD=!POSTGRES_PASSWORD!
    echo POSTGRES_DB=!POSTGRES_DB!
    echo.
    echo ENABLE_AUTH=!ENABLE_AUTH!
    echo DEFAULT_AUTH_METHOD=!DEFAULT_AUTH_METHOD_ENV!
    if "!ENABLE_AUTH!"=="1" (
        echo FLASK_SECRET_KEY=!FLASK_SECRET_KEY!
        echo LOGIN_REDIRECT_URL=!LOGIN_REDIRECT_URL!
        if /I "!DEFAULT_AUTH_METHOD!"=="auth0" (
            echo AUTH0_DOMAIN=!AUTH0_DOMAIN!
            echo AUTH0_CLIENT_ID=!AUTH0_CLIENT_ID!
            echo AUTH0_CLIENT_SECRET=!AUTH0_CLIENT_SECRET!
            echo AUTH0_CALLBACK_URL=!AUTH0_CALLBACK_URL!
            echo AUTH0_AUDIENCE=!AUTH0_AUDIENCE!
            echo AUTH0_DB_CONNECTION=!AUTH0_DB_CONNECTION!
            echo AUTH0_PUBLIC_KEY_URI=!AUTH0_PUBLIC_KEY_URI!
        )
        if /I "!DEFAULT_AUTH_METHOD!"=="shibboleth" (
            echo SHIBBOLETH_LOGIN_URL=!SHIBBOLETH_LOGIN_URL!
            echo SHIBBOLETH_LOGOUT_URL=!SHIBBOLETH_LOGOUT_URL!
        )
    )
    echo ENABLE_CHAT=0
) > %ENV_FILE%

REM Function to download and run docker-compose
:run_docker_compose
echo Setting DOCKER_COMPOSE_URL environment variable...
set DOCKER_COMPOSE_URL=https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/main/docker-compose.yml
echo Downloading docker-compose.yml from %DOCKER_COMPOSE_URL%...
curl -fsSL -o docker-compose.yml %DOCKER_COMPOSE_URL%
if errorlevel 1 (
    echo Error: Failed to download docker-compose.yml. Check the URL and try again.
    exit /b 1
)

powershell -NoProfile -Command "(Get-Content docker-compose.yml) -replace '\"\d+:5055\"','\"%APP_PORT%:5055\"' | Set-Content docker-compose.yml"
if errorlevel 1 (
    echo Error: Failed to set app port in docker-compose.yml.
    exit /b 1
)

set "MDV_MOUNT=mdv-data:/app/mdv"
if /I "%STORAGE_MODE%"=="host" set "MDV_MOUNT=%HOST_DATA_PATH%:/app/mdv"
powershell -NoProfile -Command "$mount = ($env:MDV_MOUNT -replace '\\','/'); $content = Get-Content docker-compose.yml; $content = $content -replace '^\s*-\s*.+:/app/mdv\s*$', ('      - ''{0}''' -f $mount); Set-Content docker-compose.yml $content"
if errorlevel 1 (
    echo Error: Failed to set data mount in docker-compose.yml.
    exit /b 1
)

echo Creating deployment override for port %APP_PORT%...
> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo services:
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo   mdv_app:
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo     env_file:
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo       - "%ENV_FILE%"
if /I "%POSTGRES_MODE%"=="reuse" if not "!REUSE_DB_NETWORK!"=="" (
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo     networks:
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo       - reused_db_net
)
if /I "%DB_BACKEND%"=="postgres" if /I not "%POSTGRES_MODE%"=="reuse" (
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo   mdv_db:
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo     container_name: mdv-db
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo     env_file:
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo       - "%ENV_FILE%"
)
if /I "%POSTGRES_MODE%"=="reuse" if not "!REUSE_DB_NETWORK!"=="" (
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo networks:
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo   reused_db_net:
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo     external: true
>> .docker-deploy.%DEPLOYMENT_NAME%.override.yml echo     name: !REUSE_DB_NETWORK!
)

set USE_LOCAL_MDV_IMAGE=0
docker image inspect mdvadmin/mdv:stable >nul 2>&1
if not errorlevel 1 (
    set USE_LOCAL_MDV_IMAGE=1
    echo Detected local image mdvadmin/mdv:stable. Skipping mdv_app pull.
)

if /I "%DEPLOY_MODE%"=="replace" (
    echo Stopping existing deployment "%DEPLOYMENT_NAME%"...
    docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml down --remove-orphans
    if "!RESET_DATA_VOLUME!"=="1" (
        docker volume rm %DEPLOYMENT_NAME%_mdv-data >nul 2>&1
    )
    if "!REMOVE_OLD_VOLUME!"=="1" (
        docker volume rm %DEPLOYMENT_NAME%_mdv-data >nul 2>&1
    )
    if "!REMOVE_OLD_HOST_PATH!"=="1" if not "!PREVIOUS_HOST_DATA_PATH!"=="" (
        rd /s /q "!PREVIOUS_HOST_DATA_PATH!" >nul 2>&1
    )
)

if /I "%DB_BACKEND%"=="sqlite" (
    if not exist docker-sqlite.override.yml (
        echo Error: Missing docker-sqlite.override.yml in project root.
        exit /b 1
    )

    if "!USE_LOCAL_MDV_IMAGE!"=="0" (
        echo Pulling mdv_app image...
        docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml -f .docker-deploy.%DEPLOYMENT_NAME%.override.yml -f docker-sqlite.override.yml pull mdv_app
        if errorlevel 1 (
            echo Error: Failed to pull mdv_app image.
            exit /b 1
        )
    )

    echo Running docker-compose ^(sqlite backend^)...
    docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml -f .docker-deploy.%DEPLOYMENT_NAME%.override.yml -f docker-sqlite.override.yml up -d --no-deps mdv_app
    if errorlevel 1 (
        echo Error: Failed to run docker-compose for sqlite backend.
        exit /b 1
    )

    REM Ensure postgres service is not left running from previous deploys.
    docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml stop mdv_db >nul 2>&1
    docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml rm -f mdv_db >nul 2>&1
    docker volume rm %DEPLOYMENT_NAME%_postgres-data >nul 2>&1
    if /I "!PREV_DB_BACKEND!"=="postgres" if /I "!PREV_DB_HOST!"=="mdv-db" if not "!PREV_DB_SCHEMA!"=="" (
        if "!PREV_DB_USER!"=="" set PREV_DB_USER=testuser
        if "!PREV_DB_NAME!"=="" set PREV_DB_NAME=mdv
        echo Cleaning previous postgres schema "!PREV_DB_SCHEMA!" from mdv-db...
        docker exec mdv-db psql -U "!PREV_DB_USER!" -d "!PREV_DB_NAME!" -c "DROP SCHEMA IF EXISTS ""!PREV_DB_SCHEMA!"" CASCADE;" >nul 2>&1
        set NON_SYSTEM_SCHEMA_COUNT=
        docker exec mdv-db psql -U "!PREV_DB_USER!" -d "!PREV_DB_NAME!" -t -A -c "SELECT COUNT(*) FROM information_schema.schemata WHERE schema_name NOT IN ('pg_catalog','information_schema','public') AND schema_name NOT LIKE 'pg_toast%%' AND schema_name NOT LIKE 'pg_temp_%%';" > "%TEMP%\mdv_schema_count.tmp" 2>nul
        for /f %%s in (%TEMP%\mdv_schema_count.tmp) do set NON_SYSTEM_SCHEMA_COUNT=%%s
        del /f /q "%TEMP%\mdv_schema_count.tmp" >nul 2>&1
        if "!NON_SYSTEM_SCHEMA_COUNT!"=="0" (
            echo No non-system schemas remain in mdv-db. Stopping and removing mdv-db container.
            docker stop mdv-db >nul 2>&1
            docker rm mdv-db >nul 2>&1
        )
    )
) else (
    if /I "%POSTGRES_MODE%"=="reuse" (
        if "!USE_LOCAL_MDV_IMAGE!"=="0" (
            echo Pulling mdv_app image...
            docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml -f .docker-deploy.%DEPLOYMENT_NAME%.override.yml pull mdv_app
            if errorlevel 1 (
                echo Error: Failed to pull mdv_app image.
                exit /b 1
            )
        )

        echo Running docker-compose ^(postgres reuse mode, app only^)...
        docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml -f .docker-deploy.%DEPLOYMENT_NAME%.override.yml up -d --no-deps mdv_app
        if errorlevel 1 (
            echo Error: Failed to run docker-compose for postgres reuse mode.
            exit /b 1
        )

        docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml stop mdv_db >nul 2>&1
        docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml rm -f mdv_db >nul 2>&1
    ) else (
        if "!USE_LOCAL_MDV_IMAGE!"=="1" (
            echo Pulling postgres image only ^(mdv_app is local^)...
            docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml -f .docker-deploy.%DEPLOYMENT_NAME%.override.yml pull mdv_db
            if errorlevel 1 (
                echo Error: Failed to pull postgres image. Check the configuration and try again.
                exit /b 1
            )
        ) else (
            echo Pulling stable Docker images...
            docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml -f .docker-deploy.%DEPLOYMENT_NAME%.override.yml pull
            if errorlevel 1 (
                echo Error: Failed to pull Docker images. Check the configuration and try again.
                exit /b 1
            )
        )

        echo Running docker-compose with docker-compose.yml...
        docker compose -p %DEPLOYMENT_NAME% --env-file %ENV_FILE% -f docker-compose.yml -f .docker-deploy.%DEPLOYMENT_NAME%.override.yml up -d
        if errorlevel 1 (
            echo Error: Failed to run docker-compose. Check the configuration and try again.
            exit /b 1
        )
    )
)

echo MDV application deployment completed successfully!

echo.
echo ******  Open your web browser and go to https://localhost:%APP_PORT% to access the MDV application  ******
echo.

exit /b 0
