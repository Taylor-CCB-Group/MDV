# Build the frontend with npm
# warning - we had an obscure npm error with cross-env not found, and it seems like using an earlier nodejs version fixed it
# for some reason node_modules/.bin wasn't populated after `RUN npm install` - but manually running it in the container worked
# not sure if there's a way of pinning a more specific build of the base image. May want to see if we can reproduce this issue outside of docker.
FROM nikolaik/python-nodejs:python3.12-nodejs20 AS frontend-builder

# ARG and ENV for build information, to be passed in from the CI pipeline
ARG GIT_COMMIT_DATE="unknown"
ARG GIT_BRANCH_NAME="unknown"
ARG GIT_COMMIT_HASH="unknown"
ARG GIT_LAST_COMMIT_MESSAGE="unknown"
ARG BUILD_DATE="unknown"
ARG GIT_DIRTY="false"

# Debug output to verify build args
RUN echo "Build args received:" && \
    echo "GIT_COMMIT_DATE: $GIT_COMMIT_DATE" && \
    echo "GIT_BRANCH_NAME: $GIT_BRANCH_NAME" && \
    echo "GIT_COMMIT_HASH: $GIT_COMMIT_HASH" && \
    echo "GIT_LAST_COMMIT_MESSAGE: $GIT_LAST_COMMIT_MESSAGE" && \
    echo "BUILD_DATE: $BUILD_DATE" && \
    echo "GIT_DIRTY: $GIT_DIRTY"

ENV VITE_GIT_COMMIT_DATE=$GIT_COMMIT_DATE
ENV VITE_GIT_BRANCH_NAME="$GIT_BRANCH_NAME"
ENV VITE_GIT_COMMIT_HASH=$GIT_COMMIT_HASH
ENV VITE_GIT_LAST_COMMIT_MESSAGE="$GIT_LAST_COMMIT_MESSAGE"
ENV VITE_BUILD_DATE=$BUILD_DATE
ENV VITE_GIT_DIRTY=$GIT_DIRTY

# Debug output to verify environment variables
RUN echo "Environment variables set:" && \
    echo "VITE_GIT_COMMIT_DATE: $VITE_GIT_COMMIT_DATE" && \
    echo "VITE_GIT_BRANCH_NAME: $VITE_GIT_BRANCH_NAME" && \
    echo "VITE_GIT_COMMIT_HASH: $VITE_GIT_COMMIT_HASH" && \
    echo "VITE_GIT_LAST_COMMIT_MESSAGE: $VITE_GIT_LAST_COMMIT_MESSAGE" && \
    echo "VITE_BUILD_DATE: $VITE_BUILD_DATE" && \
    echo "VITE_GIT_DIRTY: $VITE_GIT_DIRTY"

# Install system packages as root (these require root privileges)
# this layer will change less frequently than the others, so it's good to have it first
# Install HDF5 library, for some reason poetry can't install it in this context as of now
# see https://github.com/h5py/h5py/issues/2146 for similar-ish issue
RUN apt-get update && apt-get install -y libhdf5-dev || (cat /var/log/apt/term.log && exit 1)
RUN apt-get install -y netcat-openbsd || (cat /var/log/apt/term.log && exit 1)
RUN apt-get install -y telnet || (cat /var/log/apt/term.log && exit 1)
RUN apt-get install -y iputils-ping || (cat /var/log/apt/term.log && exit 1)

# Switch to the pn user early in the process to avoid permission issues
USER pn
WORKDIR /app

# Configure poetry for the pn user - let it create virtual environments for isolation
# 
# VIRTUAL ENVIRONMENTS IN DOCKER CONSIDERATIONS:
# 
# PROS of using virtual environments (current approach):
# - Security: poetry install runs in isolated environment, can't affect system packages
# - Permission isolation: avoids conflicts with system-installed packages like certifi
# - Standard practice: how poetry is intended to be used
# - Safer with untrusted code: LLM-generated code runs in more isolated environment
# 
# CONS of using virtual environments:
# - Less convenient for terminal use: need to run `poetry shell` or `poetry run` for commands
# - Seems redundant: container already provides isolation from host system
# - Extra layer: adds complexity when debugging dependency issues
# - Performance: slight overhead from virtual environment activation
# 
# PROS of disabling virtual environments (previous approach):
# - Convenience: can run python/pip commands directly in terminal without activation
# - Simpler: no need to remember to use `poetry run` or `poetry shell`
# - Direct access: easier to inspect installed packages and run ad-hoc scripts
# 
# CONS of disabling virtual environments:
# - Permission conflicts: user can't modify system packages (the current issue)
# - Security risk: poetry install can affect system Python environment
# - Mixing concerns: system packages and project packages in same namespace
# 
# DECISION: Using virtual environments for security and permission isolation,
# `.bashrc` is used to activate the virtual environment automatically.
RUN poetry config virtualenvs.in-project true

# Install Python dependencies using Poetry
# this should be early in the process because it's less likely to change
# copy poetry files first to cache the install step (with correct ownership)
COPY --chown=pn:pn python/pyproject.toml ./python/
COPY --chown=pn:pn python/poetry.lock ./python/
WORKDIR /app/python
# trouble with this is that it doesn't have the mdvtools code yet...
# still worth installing heavy dependencies here
# BUT - 
#16 20.81 Installing the current project: mdvtools (0.0.1)
#16 20.81 
#16 20.81 Warning: The current project could not be installed: [Errno 2] No such file or directory: '/app/python/README.md'
#16 20.81 If you do not want to install the current project use --no-root.
#16 20.81 If you want to use Poetry only for dependency management but not for packaging, you can disable package mode by setting package-mode = false in your pyproject.toml file.
#16 20.81 In a future version of Poetry this warning will become an error!
# seems to be ok to first install with --no-root for dependencies, then install the root package later
# we don't want to set package-mode = false in pyproject.toml
RUN poetry install --with dev,backend,auth --no-root

WORKDIR /app
# copy the package.json and package-lock.json as a separate step so npm install can be cached (with correct ownership)
COPY --chown=pn:pn package*.json ./
## Install npm dependencies
RUN npm install

# Copy only frontend-related files needed for the build
# This ensures the frontend build cache only invalidates when frontend files change
COPY --chown=pn:pn src/ ./src/
COPY --chown=pn:pn vite.config.mts ./
COPY --chown=pn:pn tsconfig.json ./
COPY --chown=pn:pn index.html ./
COPY --chown=pn:pn public/ ./public/
COPY --chown=pn:pn tailwind.config.js ./
COPY --chown=pn:pn postcss.config.js ./
COPY --chown=pn:pn biome.jsonc ./

RUN mkdir -p dist
## Run the npm build script for Flask and Vite
# This step will be cached unless frontend-related files change
RUN npm run build-flask-dockerjs

# Copy the complete project again to ensure all directories exist in final container
COPY --chown=pn:pn . .

WORKDIR /app/python
# installing again so we have mdvtools as a module, on top of the previous install layer with dependencies
# this step should be very fast
# if we don't have this, the server itself runs, but anything that doesn't run from this workdir will fail to import mdvtools
ENV NUMCODECS_NO_ASM=1
RUN poetry install --with dev,backend,auth 

# Expose the port that Flask will run on
EXPOSE 5055 

# something changed causing npm to need this in order for `source` to work in npm scripts
RUN npm config set script-shell "/bin/bash"

# Add automatic poetry shell activation for terminal use
RUN echo 'cd /app/python && $(poetry env activate)' >> ~/.bashrc && \
    echo 'cd /app' >> ~/.bashrc

#USER root
#RUN mkdir -p /app/logs && chown -R pn:pn /app/logs
#USER pn

USER root
RUN mkdir -p /app/mdv && chown -R pn:pn /app/mdv
USER pn


# Command to run Gunicorn
CMD ["poetry", "run", "gunicorn", "-k", "gevent", "-t", "0", "-w", "4", "-b", "0.0.0.0:5055", "--reload", "--capture-output", "--log-level", "info", "mdvtools.dbutils.safe_mdv_app:app"]
#CMD ["poetry", "run", "python", "-m", "mdvtools.dbutils.mdv_server_app"]


