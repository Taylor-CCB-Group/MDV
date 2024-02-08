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
FROM python:3.10.9 AS python-builder

# Set the working directory inside the container
WORKDIR /app

# Copy everything from inside the frontend-builder /app to /app
COPY --from=frontend-builder /app/python /app/python

# Set the working directory to the Python directory
WORKDIR /app/python

# Install Python dependencies
RUN pip install -e /app/python

# Set the working directory back to /app
WORKDIR /app

# Expose the port that Flask will run on
EXPOSE 5052

# Run your Python script
CMD ["python", "-m", "mdvtools.test_projects.scanpy_pbmc3k"]
