# Build the frontend with npm
# warning - we had an obscure npm error with cross-env not found, and it seems like using an earlier nodejs version fixed it
# for some reason node_modules/.bin wasn't populated after `RUN npm install` - but manually running it in the container worked
# not sure if there's a way of pinning a more specific build of the base image. May want to see if we can reproduce this issue outside of docker.
FROM nikolaik/python-nodejs:python3.12-nodejs20 AS frontend-builder

# this layer will change less frequently than the others, so it's good to have it first
# Install HDF5 library, for some reason poetry can't install it in this context as of now
# see https://github.com/h5py/h5py/issues/2146 for similar-ish issue
RUN apt-get update && apt-get install -y libhdf5-dev || (cat /var/log/apt/term.log && exit 1)
RUN apt-get install -y netcat-openbsd || (cat /var/log/apt/term.log && exit 1)
RUN apt-get install -y telnet || (cat /var/log/apt/term.log && exit 1)
RUN apt-get install -y iputils-ping || (cat /var/log/apt/term.log && exit 1)

# Install Python dependencies using Poetry
# this should be early in the process because it's less likely to change
WORKDIR /app
# copy poetry files first to cache the install step
COPY python/pyproject.toml ./python/
COPY python/poetry.lock ./python/
WORKDIR /app/python
RUN poetry config virtualenvs.create false
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
# copy the package.json and package-lock.json as a separate step so npm install can be cached
COPY package*.json ./
## Install npm dependencies
RUN npm install


# Copy the entire project to the working directory
COPY . .


# & start dev server (should only happen in dev - but other yml configs won't expose the port)
# RUN npm run dev # run manually in container

## Run the npm build script for Flask and Vite
# this will often change, so it's good to have it last... doesn't seem to be cached
RUN npm run build-flask-dockerjs


WORKDIR /app/python
# installing again so we have mdvtools as a module, on top of the previous install layer with dependencies
# this step should be very fast
# if we don't have this, the server itself runs, but anything that doesn't run from this workdir will fail to import mdvtools
RUN poetry install --with dev,backend,auth 

# Expose the port that Flask will run on
EXPOSE 5055 

# something changed causing npm to need this in order for `source` to work in npm scripts
RUN npm config set script-shell "/bin/bash"

# Command to run Gunicorn
# nb -t 0 means no timeout, this is a dev setting, chat requests are long running & should be made async
# for SocketIO, we also need to run with --worker-class eventlet or gevent https://flask-socketio.readthedocs.io/en/latest/deployment.html
# adding `-k gevent` for this
CMD ["poetry", "run", "gunicorn", "-k", "gevent", "-t", "0", "-w", "1", "-b", "0.0.0.0:5055", "--reload", "--access-logfile", "/app/logs/access.log", "--error-logfile", "/app/logs/error.log", "--capture-output", "mdvtools.dbutils.mdv_server_app:app"]
#CMD ["poetry", "run", "python", "-m", "mdvtools.dbutils.mdv_server_app"]


