
# **Simple Guide to Creating Views in MDV**

This step-by-step guide will walk you through creating an MDV (Multi-Dimensional Viewer) project, adding a table view, and then a scatterplot view. This tutorial will take you through the various steps and the full code is shown at the end. It is also available to run from the GitHub repo:
```
python python/mdvtools/test_projects/scanpy_pbmc3k_tutorial.py
```

---

## **Step 1: Setting Up the Project**

The first step is to create an MDV project from the PBMC3k dataset provided by Scanpy.

### **Code**:

```python
from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject
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

# Enable editing for the project
p.set_editable(True)
```

### **Explanation**:
1. **Set Up Project Directory**: 
   - The project is stored in `~/mdv/pbmc3k`.
   - If the directory doesnt exist, it will be created.

2. **Load or Create the Project**:
   - The PBMC3k dataset is processed and converted into an MDV project.
   - If the project already exists, it will be reused.

3. **Enable Editing**:
   - This allows modifications to the project configuration.

---

## **Step 2: Adding a Table View**

Now that the project is set up, let's add a table view to display metadata about the cells.

### **Code**:

```python
# Add table
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

    # Configure the project view with the table
    view_config = {
        "initialCharts": {
            "cells": [
                table_plot.plot_data
            ]
        }
    }
    p.set_view("default", view_config)
```

### **Explanation**:
1. **Retrieve the Data**:
   - The `cells` datasource contains information about individual cells in the dataset.

2. **Create a TablePlot**:
   - Displays metadata in a tabular format.
   - Configured to show all columns in the `cells` dataframe.

3. **Set the View**:
   - The table is added to the `cells` view in the project configuration.

---

## **Step 3: Adding a Scatterplot View**

In addition to the table, add a scatterplot view for visualising UMAP data.

### **Code**:

```python
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
```

### **Explanation**:
1. **Create a ScatterPlot**:
   - Displays a 2D scatterplot using UMAP dimensions.
   - Configured with default colors and positions.

2. **Add to View Configuration**:
   - Both the table and scatterplot are included in the `cells` view.

---

## **Step 4: Serve the Project**

Finally, serve the project to interact with the views in a web browser.

### **Code**:

```python
# Serve the project
p.serve(port=5052)
```

### **Explanation**:
- The server runs locally on port 5052.
- Open `http://localhost:5052` in your browser to view the project.

---

## **Complete Code**

```python
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
```

---

## **Conclusion**

Congratulations! You have successfully created an MDV project with a table and scatterplot view. Open your browser at `http://localhost:5052` to interact with your visualisations.

### **Next Steps**
- Experiment with adding more visualisations (e.g., histograms, 3D scatterplots).
- Link additional datasources for enhanced interactivity.
- Customise the configuration to suit your analysis needs.
