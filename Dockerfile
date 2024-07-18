# Build the frontend with npm
# warning - we had an obscure npm error with cross-env not found, and it seems like using an earlier nodejs version fixed it
# for some reason node_modules/.bin wasn't populated after `RUN npm install` - but manually running it in the container worked
# not sure if there's a way of pinning a more specific build of the base image. May want to see if we can reproduce this issue outside of docker.
FROM nikolaik/python-nodejs:python3.12-nodejs20 AS frontend-builder

# Set the working directory inside the container
WORKDIR /app

# Copy the entire project to the working directory
COPY . .

# Install npm dependencies
RUN npm install

# Run the npm build script for Flask and Vite
RUN npm run build-flask-vite

# bootstrap project folder - this won't be necessary in future
#RUN mkdir -p /app/mdv/pbmc3k /app/mdv/pbmc3k_project2

# Install HDF5 library, for some reason poetry can't install it in this context as of now
# see https://github.com/h5py/h5py/issues/2146 for similar-ish issue
RUN apt-get update && apt-get install -y libhdf5-dev

#RUN pip install gunicorn
# Install Python dependencies using Poetry
WORKDIR /app/python
RUN poetry install --with dev,backend 

# Expose the port that Flask will run on
EXPOSE 5055 

# Command to run Gunicorn
CMD ["poetry", "run", "gunicorn", "-w", "1", "-b", "0.0.0.0:5055", "mdvtools.dbutils.mdv_server_app:app"]
#CMD ["poetry", "run", "python", "-m", "mdvtools.dbutils.mdv_server_app"]


