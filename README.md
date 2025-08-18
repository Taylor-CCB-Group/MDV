# Multi-Dimensional Viewer

![logo](images/mdv_logo.png)

The Multi-Dimensional Viewer (MDV) allows visualisation and analysis of a wide range of biological data. It supports a broad range of data types, including:
•	Genomic data, such as structural variants
•	Transcriptomic data, including single-cell RNA sequencing (scRNA-Seq)
•	Proteomic data, including spatial proteomics from platforms like Hyperion and CellDive
•	Spatial transcriptomics (currently under development)

## Key Features

* Large assortment of interactive charts/widgets may be used and embedded in an analysis. These include:
    * Scatter plots (2D and 3D)
    * Box plots
    * Heat maps
    * Histograms
    * Pie charts
    * Violin plots
    * Annotation tools
    * Interactive spatial biology charts 

* Multiple views or pages of data may be created to tell a story with the data

* Charts/Widgets can pop out into separate windows to take advantage of multiple screens

* Multiple data sources (tables) can be loaded and links defined between them

* Data can be added and/or modified by the user

* Diverse range of data sources (API calls, static files) can be used by implementing custom data loaders 


* Runs in a web browser (installation not required for uploading and viewing data)


### System Requirements

- **Browser**: A modern browser (e.g., Chrome, Firefox, Safari, or Edge)
- **Memory**: At least 4 GB of RAM, even for handling large datasets (~10 million items), due to the lazy loading of data in raw bytes. When using Docker (see below) this requirement goes up to 16 or 32GB.

You can browse legacy projects at the original [MDV Website](https://mdv.molbiol.ox.ac.uk/) but we recommend looking at our new [MDV website](https://mdv.ndm.ox.ac.uk/) for resources and documentation.

### Steps to install MDV
- Install Docker from https://www.docker.com/ (Version 20.10 or Later). Docker allows the full MDV application to be installable on your local machine or server on Mac, Windows WSL2 or Linux.
- Create a folder and download the script called ‘deploy.sh’ from the Git repo:
  ```
  mkdir deploy_folder && cd deploy_folder
  wget --no-cache https://raw.githubusercontent.com/Taylor-CCB-Group/MDV/dev/deploy.sh
  ```
- Make the script executable:
  ```
  chmod u+x deploy.sh
  ```
- Run the automation script:
  ```
  ./deploy.sh
  ```
  The script will prompt to set up environment variables; provide your own values for the configuration parameters. 
  
  If you need authentication, enter 'y' when prompted. (Note: Disable it by entering 'n' if you don’t need the authentication feature.)
  
  - Once the process is completed, make sure containers are up and running using command:
    ```
    docker ps
    ```
  - At the end of installation, you will be prompted to visit http://localhost:5055 to access MDV.


You can download an example MDV project to try at:

```
https://mdv.ndm.ox.ac.uk/MDV/pbmc3k-mdv.zip
```
Then "Import an existing project" via the Project Catalogue View  
![image](https://private-user-images.githubusercontent.com/10518908/479036779-0253d7cd-6501-4f74-b6e6-27904cb65446.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3NTU1MjA2MzksIm5iZiI6MTc1NTUyMDMzOSwicGF0aCI6Ii8xMDUxODkwOC80NzkwMzY3NzktMDI1M2Q3Y2QtNjUwMS00Zjc0LWI2ZTYtMjc5MDRjYjY1NDQ2LnBuZz9YLUFtei1BbGdvcml0aG09QVdTNC1ITUFDLVNIQTI1NiZYLUFtei1DcmVkZW50aWFsPUFLSUFWQ09EWUxTQTUzUFFLNFpBJTJGMjAyNTA4MTglMkZ1cy1lYXN0LTElMkZzMyUyRmF3czRfcmVxdWVzdCZYLUFtei1EYXRlPTIwMjUwODE4VDEyMzIxOVomWC1BbXotRXhwaXJlcz0zMDAmWC1BbXotU2lnbmF0dXJlPTc2ZWE1MTVkMGE0ODZlNGIyNzk1MzY4ZTcwYzRjZTA1NzdiODY2NTZhZmU0MTBlNDNjODA3YzAzZTY3ZjQyZjMmWC1BbXotU2lnbmVkSGVhZGVycz1ob3N0In0.wnMvLr1NJCr7V1jnC_RK6712PPgoayIvyM8j7BeAZaU)

and you should see the following once uploaded.

![image](https://github.com/user-attachments/assets/e394f044-eb8a-44d3-9c05-a2576b410c19)

You can now interact with these controls or add Charts (scatterplots, rowcharts etc) and Views (contains the charts) to allow visualisation and querying of the data. A further tutorial on creating simple charts and tables is shown in [this tutorial](docs/tutorials/scanpy_python_tutorial_1.md).

See Tutorials and Documentation for how to use MDV.

## Detailed Platform Specific Installation

### Introduction

 This section guides you through configuring MDV containers on localhost with Linux/Mac.

### Prerequisites

Presuming you are a sudo user, follow these steps to ensure your system is ready:

1. **Install Docker (Version 20.10 or Later)**
   - Docker version 20.10 or later includes Docker Compose by default, so no separate installation is needed.
   - Check if Docker is installed:
     ```
     docker --version
     ```
     If installed, ensure it’s version 20.10 or later.
   - If Docker is not installed, follow these steps:
     - **Linux:**
       - For Debian/Ubuntu/CentOS/Fedora, install Docker using the official script:
         ```
         curl -fsSL https://get.docker.com | sh
         ```
       - For Arch Linux, use pacman to install Docker:
         ```
         sudo pacman -S docker && sudo systemctl enable --now docker
         ```
     - **macOS:**
       - Install Docker Desktop from the official website or via Homebrew:
         ```
         brew install --cask docker
         ```
       - If you don’t have Homebrew, install it first:
         ```
         /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
         ```

2. **Start Docker Before Running the Script**
   - Ensure Docker is running before executing the script.
   - **Linux:**
     ```
     sudo systemctl start docker
     ```
     If systemctl is unavailable, try:
     ```
     sudo dockerd
     ```
   - **macOS:**
     ```
     open -a Docker
     ```
     Or if installed via Homebrew:
     ```
     brew services start docker
     ```

3. **Ensure Ports 5055 and 5432 Are Available**
   - These ports must be free before running the script.
   - Check if the ports are in use:
     ```
     sudo lsof -i :5055
     sudo lsof -i :5432
     ```
   - Stop the process using the port:
     ```
     sudo kill -9 <PID>
     ```
     Replace `<PID>` with the actual PID from the output.

4. **Install wget**
   - Once Docker is installed and running, and the required ports are available, you’re ready to proceed!
   - **Linux (For Ubuntu/Debian-based distros):**
     ```
     sudo apt update && sudo apt install wget -y
     ```
   - **macOS:**
     Using Homebrew (Recommended):
     ```
     brew install wget
     ```


### Editing code within the Docker container

Once your MDV Docker container is up and running, you can edit the code directly within the container. Follow these steps to access and modify the code:

1. **Access the Running Container**
   - Use the following command to access the shell of the running container:
     ```
     docker exec -it <container_name> /bin/bash
     ```
     Replace `<container_name>` with the actual name or ID of your running container. You can find this by running `docker ps`.

2. **Navigate to the Code Directory**
   - Once inside the container, navigate to the directory where the MDV code is located. This is typically in a directory like `/app` or `/src`, depending on how the Dockerfile is set up.
     ```
     cd /path/to/code
     ```

3. **Edit the Code**
   - You can use command-line text editors like `vi`, `nano`, or `emacs` to edit the code. For example, to edit a file using `nano`, you would run:
     ```
     nano filename.py
     ```

4. **Save Changes and Exit**
   - After making your changes, save the file and exit the editor. For `nano`, you can do this by pressing `CTRL + O` to save and `CTRL + X` to exit.

5. **Restart the Container (if necessary)**
   - If your changes require a restart of the application, you can restart the container using:
     ```
     docker restart <container_name>
     ```

This process allows you to make and test changes in real time within the Docker environment. Ensure you have the necessary permissions to edit files within the container.

### Editing code with VSCode and Cursor

In addition to command-line editors, you can use modern IDEs like VSCode and Cursor to edit code within the Docker container. Here’s how:

#### Using VSCode (Dev Containers)

1. **Install the Dev Containers extension**
   - Install the [Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension in VSCode.

2. **Open the project in a Dev Container**
   - With the repository open in VSCode, use the command palette (`Ctrl+Shift+P` or `Cmd+Shift+P` on Mac) and select `Dev Containers: Reopen in Container`.
   - If prompted, VSCode will build the container using the project's configuration and reopen the workspace inside the container.

3. **Edit the code**
   - Once attached, you can navigate and edit the code as you would in a local environment. All VSCode features, such as IntelliSense and debugging, are available.


These options provide a more graphical and feature-rich environment for editing code within Docker containers, enhancing productivity and ease of use.
