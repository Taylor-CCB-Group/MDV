from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject
import scanpy as sc
import os

base = os.path.expanduser('~/mdv')
project_folder = os.path.join(base, 'pbmc3k')
if not os.path.exists(os.path.expanduser('~/mdv')):
    os.makedirs(base)
if not os.path.exists(project_folder):
    data = sc.datasets.pbmc3k_processed()
    p = convert_scanpy_to_mdv(project_folder, data)
else:
    print('using existing project...')
    p = MDVProject(project_folder)
p.set_editable(True)
p.serve(port=5052) # port conflict locally as of writing...