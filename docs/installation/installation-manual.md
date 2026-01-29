# MDV Installation on Linux, Mac, Windows (WSL2)

Welcome to the MDV setup manual.

This document guides you through configuring MDV containers, on localhost with Linux, Mac or Windows (WSL2).

## Prerequisites

Presuming you are a sudo user, follow these steps to ensure your system is ready:

### 1. Install Docker (Version 20.10 or Later)

Docker v20.10+ includes Docker Compose by default, so no separate installation is required.

**Check if Docker is installed:**

Run the following command (works on Linux, macOS, and Windows):

```bash
docker --version
```

You should see output similar to:

```
Docker version 20.10.24, build 297e128
```

- If installed, ensure it's version 20.10 or later.

**If Docker is not installed, follow these steps:**

**Windows 10/11 (Pro, Enterprise, or Home with WSL2):**

- Download [Docker Desktop for Windows](https://www.docker.com/products/docker-desktop/)
- Install and restart your system.
- Ensure WSL2 (Windows Subsystem for Linux) and Virtualization are enabled:
  ```bash
  wsl --install
  ```
- Verify Docker Desktop is running in the system tray.

**Linux:**

For Debian/Ubuntu/CentOS/Fedora:

```bash
curl -fsSL https://get.docker.com | sh
```

For Arch Linux:

```bash
sudo pacman -S docker && sudo systemctl enable --now docker
```

**macOS (Intel & Apple Silicon):**

- Install Docker Desktop from the [official website](https://www.docker.com/products/docker-desktop/)
- Or install via Homebrew:
  ```bash
  brew install --cask docker
  ```
- If you don't have Homebrew, install it first:
  ```bash
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  ```

### 2. Start Docker Before Running the Script

Docker must be running before executing the MDV setup script.

**Windows (Docker Desktop with WSL2):**

- Ensure Docker Desktop is running (check system tray).
- If using WSL, start Ubuntu with:
  ```bash
  wsl -d Ubuntu
  ```
  (This ensures the WSL2 backend is active.)

**macOS (Docker Desktop):**

- Start Docker from Applications folder, or run:
  ```bash
  open -a Docker
  ```
- If installed via Homebrew:
  ```bash
  brew services start docker
  ```

**Linux:**

```bash
sudo systemctl start docker
```

If systemctl is unavailable:

```bash
sudo dockerd
```

### 3. Ensure Ports 5055 and 5432 Are Available

MDV requires port 5055 (application) and port 5432 (PostgreSQL).

Before running the deployment script, make sure these ports are free.

**Check if the ports are in use:**

**Windows (PowerShell / WSL):**

Check if the ports are in use:

```bash
netstat -ano | findstr :5055
netstat -ano | findstr :5432
```

If a process is using these ports, you'll see an output like this:

```
COMMAND   PID  USER   FD   TYPE DEVICE SIZE/OFF NODE NAME 
python   1234 user   10u  IPv4  12345      0t0  TCP *:5055 (LISTEN)
```

Terminate process if occupied:

```bash
taskkill /PID <PID> /F
```

If using WSL2 (Ubuntu) backend:

```bash
sudo lsof -i :5055
sudo lsof -i :5432
sudo kill -9 <PID>
```

**macOS:**

Check if the ports are in use:

```bash
sudo lsof -i :5055
sudo lsof -i :5432
```

Terminate process if occupied:

```bash
sudo kill -9 <PID>
```

**Linux:**

Check if the ports are in use:

```bash
sudo lsof -i :5055
sudo lsof -i :5432
```

Terminate process if occupied:

```bash
sudo kill -9 <PID>
```

Alternative for older distros without lsof:

```bash
sudo netstat -tulpn | grep :5055
sudo netstat -tulpn | grep :5432
sudo kill -9 <PID>
```

### 4. Install wget

The MDV deployment script requires a command-line tool to download files (wget or curl).

**Windows:**

No wget installation required — Windows 10/11 includes curl by default.

**macOS:**

Install wget via Homebrew:

```bash
brew install wget
```

**Linux:**

Debian / Ubuntu:

```bash
sudo apt update && sudo apt install wget -y
```

Fedora / CentOS / RHEL:

```bash
sudo dnf install wget -y
```

Arch Linux / Manjaro:

```bash
sudo pacman -S wget
```

**Tips:**

Verify installation:

```bash
wget --version   # For macOS / Linux
curl --version   # For Windows
```

## Steps to Install MDV

Once the prerequisites are complete (Docker installed and running, ports free, wget or curl available), follow these steps to install MDV.

**Linux/macOS:**

If you are on Linux or macOS, directly jump to Step 1.

**Windows:**

Start WSL (Windows Subsystem for Linux)

Windows users (WSL2 Ubuntu):

Open PowerShell or Windows Terminal and start Ubuntu:

```bash
wsl -d Ubuntu
```

This opens an Ubuntu shell where all following commands should be executed, now go to Step 1.

### Step 1: Create a Deployment Folder and Download the Script

Run the following commands in your terminal (Linux/macOS/WSL Ubuntu):

```bash
mkdir deploy_folder && cd deploy_folder
```

Download the deployment script:

```bash
wget --no-cache https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/main/deploy.sh
```

**Tip:** On Windows WSL, wget must be installed inside Ubuntu. If not:

```bash
sudo apt update && sudo apt install wget -y
```

### Step 2: Make the Script Executable

```bash
chmod u+x deploy.sh
```

### Step 3: Run the Deployment Script

```bash
./deploy.sh
```

You will be prompted to configure environment variables:

- Enter new values or press Return to keep existing defaults as shown below.

**Configuration options:**

- **Authentication (auth_enable):** `y` → enable, `n` → disable
- **Chat support (chat_enable):** `y` → enable (requires OPENAI_API_KEY), `n` → disable

⚠️ **Password inputs are hidden for security.**

### Step 4: Verify Containers Are Running

```bash
docker ps
```

Ensure that MDV services are listed and running.

### Step 5: Access the MDV Portal

Open a browser and visit:

```
http://localhost:5055
```

✅ If everything is running correctly, the MDV portal should appear.

## Copying a Project from Host to Container

MDV uses the user `pn` from the base image with UID 1000 and GID 1000.

**Note:** On Windows, run all commands inside the WSL Ubuntu terminal (`wsl -d Ubuntu`).

**Steps:**

1. **Copy the project into the container:**

   ```bash
   sudo docker cp <host_project_path> <container_id>:/app/mdv/
   ```

2. **Set correct ownership inside the container:**

   First, list running containers:

   ```bash
   docker ps
   ```

   Identify the container ID of `mdvapp` (do not use the db container).

   Then run:

   ```bash
   docker exec -u 0 -it <mdvapp_container_id> chown -R 1000:1000 /app/mdv/<project_name>
   ```

3. **Rescan projects via MDV portal:**

   Open in your browser:

   ```
   http://localhost:5055/rescan_projects
   ```

   ⚠️ If `/rescan_projects` does not work, restart the container:

   ```bash
   docker restart <container_id>
   ```

