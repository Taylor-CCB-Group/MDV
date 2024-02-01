from mdvtools.conversions import convert_scanpy_to_mdv
import scanpy as sc
import os

project_folder = os.path.expanduser('~/mdv/pbmc3k')
data = sc.datasets.pbmc3k_processed()
p = convert_scanpy_to_mdv(project_folder, data)
p.set_editable(True)
p.serve(port=5052) # port conflict locally as of writing...