# Stage 1: Install Node.js
FROM node:16.13.1 AS frontend-builder

# Set the working directory inside the container
WORKDIR /app

# Copy package.json and package-lock.json to the working directory
COPY package.json .
COPY package-lock.json .

# Install npm dependencies
RUN npm install

# Copy the entire project to the working directory
COPY . .

# Stage 2: Build the Python backend
FROM python:3.12 AS python-builder

# Set the working directory inside the container
WORKDIR /app

# Copy everything from inside the frontend-builder /app to /app
COPY --from=frontend-builder /app/python /app/python

# Install Python dependencies
RUN npm run poetry-setup
RUN npm run python-setup

# Expose the port that Flask will run on
EXPOSE 5055

# Run your Python script
CMD ["python", "-m", "mdvtools.dbutils.mdv_server_app"]

