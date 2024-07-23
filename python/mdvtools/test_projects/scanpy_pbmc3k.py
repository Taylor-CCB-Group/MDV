from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject
import scanpy as sc
import os

base = os.path.expanduser("~/mdv")
project_folder = os.path.join(base, "pbmc3k")
if not os.path.exists(os.path.expanduser("~/mdv")):
    os.makedirs(base)
if not os.path.exists(project_folder):
    data = sc.datasets.pbmc3k_processed()
    p = convert_scanpy_to_mdv(project_folder, data)
else:
    print("using existing project...")
    p = MDVProject(project_folder)
p.set_editable(True)


def setup_views():
    """
    This is (will be) an example of how to build a set of views that are relevant to the given data...
    We could go into a bit of introspection on the project datasources & grouping columns based on what we can infer
    from names like "X_umap_1" "X_umap_2" being suited to making a 2D scatterplot...
    i.e. parse the 'cells' datasource columns & make charts for
    - 2D UMAP
    - 2D tsne
    - 3D pca
    All of those should be coloured by 'louvain' which we know from experience will show the clusters.
    The tooltip could be 'cell_id' so that the user can explore how particular gene sequences are distributed in those clusters.
    A table chart alongside that (perhaps not including the UMAP/tsne/PCA columns) can allow the user to see which values are
    represented as they interactively explore the data.
    This 'view' might only include 'cells'.

    There is also a 'link' in the project, which means that charts can be created dynamically in 'cells' based on interactions with 'genes'.
    Some of the UX in this area, as well as the APIs, are subject to review in the near-term... but a second 'view' could be programmatically
    added with (for example) an empty 'cells' and a 3D scatterplot of 'PCs' & table in 'genes'.
    Then the '+' icon ('Add Gene Expr') in the 'cells' toolbar can be used by the user at runtime to populate 'cells' with charts
    corresponding to selections made in 'genes'.
    """
    global p
    # this may not be the most efficient method, but we something in a familiar (and typed) format, so that's nice
    cell_df = p.get_datasource_as_dataframe("cells")
    # the rest is left as an exercise for the reader for now...
    print(cell_df)


p.serve(port=5052)  # port conflict locally as of writing...
