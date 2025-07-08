from .conversions import convert_scanpy_to_mdv
from .benchmarking.memory import log_memory_usage
import scanpy as sc
import argparse
import os

"""
Utility script to convert AnnData files to MDV format.

Reads an AnnData file in h5ad format and converts it to MDV, with sane defaults for large datasets.
"""



def get_file_size_gb(filepath):
    """Get file size in GB."""
    return os.path.getsize(filepath) / (1024**3)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert AnnData to MDV format")
    parser.add_argument("input_h5ad", help="Input h5ad file")
    parser.add_argument("output_folder", help="Output folder for MDV project")
    parser.add_argument("--delete-existing", action="store_false", default=True, help="Delete existing project data")
    # parser.add_argument("--add-layers", default=False, help="Add additional matrix layers (expensive)")
    
    args = parser.parse_args()

    file_size_gb = get_file_size_gb(args.input_h5ad)
    print(f"File size: {file_size_gb:.2f} GB")
    

    print(f"Loading data from {args.input_h5ad}")
    log_memory_usage("before loading")
    
    adata = sc.read_h5ad(args.input_h5ad, backed="r")
    print(f"Loaded data shape: {adata.shape}")
    log_memory_usage("after loading")
    
    import time
    print("Loading complete, converting to MDV")
    start_time = time.time()
    convert_scanpy_to_mdv(
        args.output_folder, 
        adata, 
        chunk_data=True,
        add_layer_data=False,
        delete_existing=args.delete_existing
    )
    elapsed = time.time() - start_time
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)
    print(f"Conversion complete (time spent: {minutes} min {seconds} sec)")
    log_memory_usage("after conversion")