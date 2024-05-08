# Stage 1: Build the frontend with npm
FROM node:14 AS frontend-builder

# Set the working directory inside the container
WORKDIR /app

# Copy package.json and package-lock.json to the working directory
COPY package.json .
COPY package-lock.json .

# Install npm dependencies
RUN npm install

# Copy the entire project to the working directory
COPY . .

# Run the npm build script for Flask and Vite
RUN npm run build-flask-vite

# Stage 2: Build the Python backend
FROM python:3.12 AS python-builder

# Set the working directory inside the container
WORKDIR /app

# Install HDF5 library
RUN apt-get update && apt-get install -y libhdf5-dev

# Copy the Python backend source code
COPY --from=frontend-builder /app/python /app/python

# Install Poetry
RUN curl -sSL https://install.python-poetry.org | python3 -

# Set up the PATH environment variable to include Poetry's bin directory
ENV PATH="${PATH}:/root/.local/bin"

# Copy pyproject.toml and poetry.lock from the source directory
COPY /python/pyproject.toml /python/poetry.lock /app/python/


# Set Poetry to use the specified Python interpreter
ENV POETRY_PYTHON="/usr/local/bin/python"

# Install Python dependencies using Poetry
WORKDIR /app/python
RUN poetry install --with dev,backend 

# Print list of installed packages using pip
RUN pip list

ENV PATH="/usr/local/bin:${PATH}"

# Expose the port that Flask will run on
EXPOSE 5055 

# Set the working directory to the Python directory
WORKDIR /app/python

# Run your Python script
CMD ["/usr/local/bin/python", "-m", "mdvtools.dbutils.mdv_server_app"]
