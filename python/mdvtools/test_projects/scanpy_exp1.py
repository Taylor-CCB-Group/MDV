from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject
import scanpy as sc
import os
import shutil
import mdvtools.chart_prototypes as cp

#base = os.path.expanduser('~/mdv')
#project_folder = os.path.join(base, 'pbmc3k')
#if not os.path.exists(os.path.expanduser('~/mdv')):
#    os.makedirs(base)
#if not os.path.exists(project_folder):
#    data = sc.datasets.pbmc3k_processed()
#    p = convert_scanpy_to_mdv(project_folder, data)
#else:
#    print('using existing project...')
#    p = MDVProject(project_folder)
#p.set_editable(True)


#p.serve(port=5052) # port conflict locally as of writing...


def copy_files_to_nfs(source_dir, nfs_server, nfs_dir, project_name):
    # Path to project folder on NFS server
    nfs_project_folder = os.path.join(nfs_dir, project_name)  
    
    # Create the project folder if it doesn't exist
    os.makedirs(nfs_project_folder, exist_ok=True)

    # Iterate over files and directories in the source directory
    for item in os.listdir(source_dir):
        source_item = os.path.join(source_dir, item)
        dest_item = os.path.join(nfs_project_folder, item)

        # If the item is a file and already exists in the destination directory
        if os.path.isfile(source_item) and os.path.exists(dest_item):
            # Generate a new filename with version number
            dest_item_v2 = generate_versioned_filename(dest_item)
            # Copy the file to the destination directory with the new name
            shutil.copyfile(source_item, dest_item_v2)
            print(f"File '{item}' already exists. Renamed as '{os.path.basename(dest_item_v2)}'.")
        # If the item is a directory, recursively copy its contents
        elif os.path.isdir(source_item):
            copy_files_to_nfs(source_item, nfs_server, nfs_dir, project_name)
            print(f"Directory '{item}' copied successfully.")
        else:
            # Copy the file to the destination directory
            shutil.copy2(source_item, dest_item)
            print(f"File '{item}' copied successfully.")

def generate_versioned_filename(filename):
    """
    Generates a new filename with a version number.
    Example: filename.txt -> filename_v2.txt
    """
    base_name, ext = os.path.splitext(filename)
    version = 1
    while os.path.exists(f"{base_name}_v{version}{ext}"):
        version += 1
    return f"{base_name}_v{version}{ext}"

# Example usage:
source_dir = '/Users/jayesh/mdv/pbmc3k'
nfs_server = 'localhost'  # Assuming NFS server is localhost
nfs_dir = '/Users/jayesh/nfs_www'  # Assuming the NFS share is mounted on this directory
project_name = 'pbmc3k_project2'

copy_files_to_nfs(source_dir, nfs_server, nfs_dir, project_name)
