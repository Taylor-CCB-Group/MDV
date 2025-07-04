import scanpy as sc
from mdvtools.conversions import convert_scanpy_to_mdv

# project_path = "mdv/544"
#dataset_path = "mdv/544/TAURUS_raw_counts_annotated_final.h5ad"

# adata = sc.read_h5ad(dataset_path)

#project = convert_scanpy_to_mdv(project_path, adata)

#project.serve()

# print(adata)

# print(adata.obsm)



import scanpy as sc
import pandas as pd

# Load the h5ad file
dataset_path = "/Users/mariak/Documents/TAURUS_raw_counts_annotated_final.h5ad"
adata = sc.read_h5ad(dataset_path)

# Load your coordinates (replace with your actual file path)
df = pd.read_csv("/Users/mariak/Documents/UMAP_combined_objects.txt", sep="\t", index_col=0, header=None)  # or header=[0] if it has column names

# Check that the row indices match adata.obs_names
# If necessary, align:
df = df.loc[df.index.isin(adata.obs_names)]

df = df.reindex(adata.obs_names)

# Optional: Check if all cells are covered
assert all(df.index == adata.obs_names), "Mismatch in cell ordering — consider reindexing"

# Add to obsm
adata.obsm["X_UMAP"] = df.values  # or .to_numpy()

# Save updated h5ad
adata.write_h5ad("/Users/mariak/Documents/TAURUS_raw_counts_annotated_final_UMAP.h5ad")





