from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.table_plot import TablePlot
import scanpy as sc
import os

# Set up project directory
base = os.path.expanduser("~/mdv")
project_folder = os.path.join(base, "pbmc3k")
if not os.path.exists(base):
    os.makedirs(base)

# Load data or create a new project
if not os.path.exists(project_folder):
    data = sc.datasets.pbmc3k_processed()
    p = convert_scanpy_to_mdv(project_folder, data)
else:
    print("Using existing project...")
    p = MDVProject(project_folder)

p.set_editable(True)

# Add visualisations
def setup_views():
    global p
    # Get the cells dataframe
    cell_df = p.get_datasource_as_dataframe("cells")

    # Add a table for displaying metadata
    table_plot = TablePlot(
        title="Metadata Table",
        params=list(cell_df.columns),
        size=[600, 500],
        position=[850, 10],
    )

    # Add a scatterplot for UMAP visualization
    umap_plot = ScatterPlot(
        title="UMAP 2D Visualisation",
        params=["X_umap_1", "X_umap_2"],
        size=[400, 400],
        position=[10, 10],
        default_color="#377eb8",
    )
    
    # Configure and add views to the project
    view_config = {
        "initialCharts": {
            "cells": [
                table_plot.plot_data,
                umap_plot.plot_data
            ]
        }
    }
    p.set_view("default", view_config)

# Call the setup_views function to add the visualizations
setup_views()

# Serve the project
p.serve(port=5052)
