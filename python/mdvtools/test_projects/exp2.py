import os
import paramiko

def copy_files_to_nfs(source_dir, nfs_server, nfs_dir, project_name, username, password):
    # Create SSH client
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    # Connect to SSH server with password authentication
    ssh.connect(nfs_server, username=username, password=password)

    # Path to project folder on NFS server
    nfs_project_folder = os.path.join(nfs_dir, project_name)

    # Create the project folder if it doesn't exist
    ssh.exec_command(f'mkdir -p {nfs_project_folder}')

    # Iterate over files and directories in the source directory
    for item in os.listdir(source_dir):
        source_item = os.path.join(source_dir, item)
        dest_item = os.path.join(nfs_project_folder, item)

        # If the item is a file and already exists in the destination directory
        if os.path.isfile(source_item) and file_exists_on_server(ssh, dest_item):
            # Generate a new filename with version number
            dest_item_v2 = generate_versioned_filename(ssh, dest_item)
            # Copy the file to the destination directory with the new name
            copy_file_to_server(source_item, ssh, dest_item_v2)
            print(f"File '{item}' already exists. Renamed as '{os.path.basename(dest_item_v2)}'.")
        # If the item is a directory, recursively copy its contents
        elif os.path.isdir(source_item):
            copy_files_to_nfs(source_item, nfs_server, nfs_dir, project_name, username, password)
            print(f"Directory '{item}' copied successfully.")
        else:
            # Copy the file to the destination directory
            copy_file_to_server(source_item, ssh, dest_item)
            print(f"File '{item}' copied successfully.")

    ssh.close()

def file_exists_on_server(ssh, filepath):
    """
    Check if a file exists on the server.
    """
    stdin, stdout, stderr = ssh.exec_command(f'test -e {filepath} && echo "1" || echo "0"')
    return stdout.read().strip() == b'1'

def copy_file_to_server(local_file, ssh, remote_file):
    """
    Copy a file from local machine to remote server.
    """
    sftp = ssh.open_sftp()
    sftp.put(local_file, remote_file)
    sftp.close()

def generate_versioned_filename(ssh, filename):
    """
    Generates a new filename with a version number on the server.
    Example: filename.txt -> filename_v2.txt
    """
    base_name, ext = os.path.splitext(filename)
    version = 1
    while file_exists_on_server(ssh, f"{base_name}_v{version}{ext}"):
        version += 1
    return f"{base_name}_v{version}{ext}"

# Example usage:
source_dir = '/Users/jayesh/mdv/pbmc3k'
nfs_server = 'localhost'  # Assuming NFS server is localhost
nfs_dir = '/Users/jayesh/nfs_www'  # Assuming the NFS share is mounted on this directory
project_name = 'pbmc3k_project2'
username = 'jayesh'  # Replace 'your_username' with your actual NFS server username
password = os.environ.get('NFS_PASSWORD')  # Retrieve password from environment variable

copy_files_to_nfs(source_dir, nfs_server, nfs_dir, project_name, username, password)
