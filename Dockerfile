# Build the frontend with npm
# warning - we had an obscure npm error with cross-env not found, and it seems like using an earlier nodejs version fixed it
# for some reason node_modules/.bin wasn't populated after `RUN npm install` - but manually running it in the container worked
# not sure if there's a way of pinning a more specific build of the base image. May want to see if we can reproduce this issue outside of docker.
FROM nikolaik/python-nodejs:python3.12-nodejs20 AS frontend-builder

# Set the working directory inside the container
WORKDIR /app

# Copy the entire project to the working directory
COPY . .

### !!! HOTFIX !!! ###
## There is an issue with the frontend build that appears to originate from something
## somewhat spooky - in that as far as we know, commits that had previously deployed
## successfully are now failing. 
## This is a hotfix to use a build output that is generated with
## `npm run build-flask-docker-hotfix` in an envrionment where the build is working.
## this is included in the repo for now, but should be removed once the issue is resolved.

## Install npm dependencies
# RUN npm install

## Run the npm build script for Flask and Vite
# RUN npm run build-flask-dockerjs

RUN ln -s /app/python/dist_hotfix /app/dist
### !!! /HOTFIX !!! ###

# bootstrap project folder - this won't be necessary in future
#RUN mkdir -p /app/mdv/pbmc3k /app/mdv/pbmc3k_project2

# Install HDF5 library, for some reason poetry can't install it in this context as of now
# see https://github.com/h5py/h5py/issues/2146 for similar-ish issue
RUN apt-get update && apt-get install -y libhdf5-dev || (cat /var/log/apt/term.log && exit 1)
RUN apt-get install -y netcat-openbsd || (cat /var/log/apt/term.log && exit 1)
RUN apt-get install -y telnet || (cat /var/log/apt/term.log && exit 1)
RUN apt-get install -y iputils-ping || (cat /var/log/apt/term.log && exit 1)

#RUN pip install gunicorn
# Install Python dependencies using Poetry
WORKDIR /app/python
RUN poetry install --with dev,backend 

# Expose the port that Flask will run on
EXPOSE 5055 

# Command to run Gunicorn
CMD ["poetry", "run", "gunicorn", "-w", "1", "-b", "0.0.0.0:5055", "--reload", "--access-logfile", "/app/logs/access.log", "--error-logfile", "/app/logs/error.log", "--capture-output", "mdvtools.dbutils.mdv_server_app:app"]
#CMD ["poetry", "run", "python", "-m", "mdvtools.dbutils.mdv_server_app"]


