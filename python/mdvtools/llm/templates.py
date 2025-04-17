from mdvtools.mdvproject import MDVProject

prompt_data = """
Your task is to:  
1. Identify the type of data the user needs (e.g., categorical, numerical, etc.).  
2. Use only the two DataFrames provided:
   - df1: cells (data_frame_obs)
   - df2: genes (data_frame_var)
3. Column selection logic:
   - For non-gene queries: select columns from df1 only. Inspect df1, using df1.columns
   - For gene-related queries (e.g., expression of a gene, comparison of genes, highest expressing genes):
       a. Use ONLY gene names from df2["name"] — do NOT use gene IDs or any other columns (e.g., df2["gene_ids"]).
       b. If a specific gene is mentioned by the user, check if it exists in df2["name"].
           - If it exists, return it.
           - If it does not exist, ignore it and select one or more genes from df2["name"].
       c. If no gene is mentioned, select one or more gene names from df2["name"].
       d. Only use values from df2["name"] — do NOT use any other columns from df2.
4. Always return the list of required columns as a quoted comma-separated string, like:
   - "col1", "col2"
   - Or for gene-related: "col", "gene_name"   (make sure "col" is from df1) 
5. For gene-related queries:
   - Return both df1 columns and the selected gene name (from df2["name"]).
   - Only return the name as a string (e.g., "gene_name")—do not wrap it.
6. NEVER create new DataFrames or modify existing ones.
7. Ensure that the selected columns match the visualization requirements:  
    - Abundance Box plot: Requires three categorical columns.  
      - If only one categorical variable is available, return it three times.  
      - If two are available, return one of them twice.  
    - Box plot: Requires one categorical column and one numerical column.  
    - Density Scatter plot: Requires two numerical columns and one categorical column.  
    - Dot plot: Requires one categorical column and any number of numerical columns.  
    - Heat map: Requires one categorical column and any number of numerical columns.  
    - Histogram: Requires one numerical column.  
    - Multiline chart: Requires one numerical column and one categorical column.  
    - Pie Chart: Requires one categorical column.  
    - Row Chart: Requires one categorical column.  
    - Row summary box: Requires any column(s).  
    - Sankey diagram: Requires two categorical columns.  
      - If only one categorical variable is available, return it twice.  
    - Scatter plot (2D): Requires two numerical columns and one any column for color.
    - Scatter plot (3D): Requires three numerical columns and one any column for color.  
    - Selection dialog plot: Requires any column.  
    - Stacked row chart: Requires two categorical columns.  
      - If only one categorical variable is available, return it twice.  
    - Table: Requires any column(s).  
    - Text box: Requires no columns, just text.  
    - Violin plot: Requires one categorical column and one numerical column.  
    - Wordcloud: Requires one categorical column.  
8. Important: Clearly separate the selected columns with quotes and commas.
9. NEVER explain your reasoning. Only return the required columns as a string: `"col1", "col2", "col3"`
"""

prompt_data_original = """
Your task is to:  
1. Identify the type of data the user needs (e.g., categorical, numerical, etc.).  
2. Select the most relevant column names from the first DataFrame provided unless handling a gene-related query.  
3. If the query is gene-related (e.g., gene expression value, most expressing gene, target expression, etc.), retrieve gene identifiers from the second DataFrame by:
  - Identifying the column that contains gene identifiers (the column named "name").
  - If the user specifies a gene identifier, check if it exists in the "name" column:
      - If it exists, use it.
      - If it does not, select a valid gene identifier from that column.
  - Always return the selected gene identifier along with the relevant columns from the first DataFrame.
4. Do NOT create new DataFrames. Always use the existing ones provided.
5. Ensure that the selected columns match the visualization requirements:  
    - Abundance Box plot: Requires three categorical columns.  
      - If only one categorical variable is available, return it three times.  
      - If two are available, return one of them twice.  
    - Box plot: Requires one categorical column and one numerical column.  
    - Density Scatter plot: Requires two numerical columns and one categorical column.  
    - Dot plot: Requires one categorical column and any number of numerical columns.  
    - Heat map: Requires one categorical column and any number of numerical columns.  
    - Histogram: Requires one numerical column.  
    - Multiline chart: Requires one numerical column and one categorical column.  
    - Pie Chart: Requires one categorical column.  
    - Row Chart: Requires one categorical column.  
    - Row summary box: Requires any column(s).  
    - Sankey diagram: Requires two categorical columns.  
      - If only one categorical variable is available, return it twice.  
    - Scatter plot (2D): Requires two numerical columns.  
    - Scatter plot (3D): Requires three numerical columns and one categorical column for color.  
    - Selection dialog plot: Requires any column.  
    - Stacked row chart: Requires two categorical columns.  
      - If only one categorical variable is available, return it twice.  
    - Table: Requires any column(s).  
    - Text box: Requires no columns, just text.  
    - Violin plot: Requires one categorical column and one numerical column.  
    - Wordcloud: Requires one categorical column.  
6. Return the column names in a string format, e.g., "col1", "col2", unless the query is gene-related.
7. If the query is gene-related, return the selected columns from the first DataFrame and the selected gene identifier 
from the second DataFrame in a string format, e.g., "col1", "col2", "gene name".
8. Do not provide additional explanations—only return the string.
9. ALWAYS return a string of columns, even when the query does not define them.
"""



packages_functions = """import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.charts.density_scatter_plot import DensityScatterPlot
from mdvtools.charts.heatmap_plot import HeatmapPlot
from mdvtools.charts.histogram_plot import HistogramPlot
from mdvtools.charts.dot_plot import DotPlot
from mdvtools.charts.box_plot import BoxPlot
from mdvtools.charts.scatter_plot_3D import ScatterPlot3D
from mdvtools.charts.row_chart import RowChart
from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.abundance_box_plot import AbundanceBoxPlot
from mdvtools.charts.stacked_row_plot import StackedRowChart
from mdvtools.charts.ring_chart import RingChart
from mdvtools.charts.pie_chart import PieChart
from mdvtools.charts.violin_plot import ViolinPlot
from mdvtools.charts.multi_line_plot import MultiLinePlot
from mdvtools.charts.table_plot import TablePlot
from mdvtools.charts.wordcloud_plot import WordcloudPlot
from mdvtools.charts.text_box_plot import TextBox
from mdvtools.charts.row_summary_box_plot import RowSummaryBox
from mdvtools.charts.selection_dialog_plot import SelectionDialogPlot
from mdvtools.charts.sankey_plot import SankeyPlot

import json
import numpy as np
import sys
"""
# def load_data(path):
#     #Load data from the specified CSV file.
#     return pd.read_csv(path, low_memory=False)

# def convert_plot_to_json(plot):
#     #Convert plot data to JSON format.
#     return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))


def get_createproject_prompt_RAG(project: MDVProject, path_to_data: str, datasouce_name: any, final_answer: str, question: str) -> str:
    """
    Returns the prompt for the create project RAG.
    
    Args:
    project: MDVProject we pass in the project object to get the project directory - previously we were passing the project id
      and making (wrong) assumptions about the project directory based on that - while also being less.
    path_to_data: str, path to an anndata object, subject to change
    datasouce_name: any, we need to better clarify how we reason about multi-ds queries etc...
    final_answer: str, output from the dataframe agent, which should contain a list of column names to be used in generating the view
    """
    prompt_RAG = (
        """
Context: {context}

The collection of Python scripts provided in the context, is designed to generate various types of data visualizations
using the mdvtools library. Each script focuses on a specific type of plot and follows a common structure that includes loading
data from a file, creating a plot using specific parameters, and serving the visualization through an MDV project.

All scripts in the context share a common workflow:

Setup: Define the project path, data path, and view name, the project path should always be: project_path = '"""
        + project.dir
        + """'
Plot function definition: Define the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) using a function in the same way as the context.
Project Creation: Initialize an MDVProject instance using the method: MDVProject(project_path, delete_existing=True).
Data Loading: Load data from the specified file into a pandas DataFrame using the load_data(data_path) function.
Data adding: Add the data source to the project using the method: project.add_datasource(datasource_name, data).
Plot Creation: Create the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) and define the plot paramaters in the same way as in the context.
Data Conversion: Convert the plot data to JSON format for integration with the MDV project using the convert_plot_to_json(plot) function.
Serving: Configure the project view, set it to editable, and serve the project using the .set_view(view_name, plot_view), .set_editable(True) and .serve() methods.

If no script is relevant, guided by the context generate a new script.

This string """
        + final_answer
        + """ specifies the names of the data fields that need to be plotted, for example in the params field. Get the structure of params definition from the context.

IMPORTANT: All string interpolations MUST use Python formatted string literals (f-strings).

Example:
# Correct:
message = f"gs|{{param2}}(gs)|{{param2_index}}"

# Common mistake to avoid:
WRONG: message = "gs|{{param2}}(gs)|{{param2_index}}"

The question for you to try to plot a graph based on the above instructions is given by this question: """ + question + """

If the question asks for a gene-related query (e.g., gene expression value, most expressing gene, target expression, etc.), make sure that:
1. You load both datasources that are needed, e.g. cells and genes.
2. If gene expression is required, make sure the gene id, is given as a param in the formatted string literal format: f"gs|{{param2}}(gs)|{{param2_index}}", with param2 and param2_index given by param2 = "param2"
and param2_index = data_frame_var.index.get_loc(param2) respecitvely.

The data_path are given by this variable `""" + path_to_data + """`
The datasource_name is given by this variable `""" + datasouce_name + """`

"""
    )
    return prompt_RAG

# The data should be loaded in the same way as in the examples.


def get_createproject_prompt_RAG_attempt1(project: MDVProject, path_to_data: str, datasouce_name: any, final_answer: str, question: str) -> str:
    """
    Returns the prompt for the create project RAG.
    
    Args:
    project: MDVProject we pass in the project object to get the project directory - previously we were passing the project id
      and making (wrong) assumptions about the project directory based on that - while also being less.
    path_to_data: str, path to an anndata object, subject to change
    datasouce_name: any, we need to better clarify how we reason about multi-ds queries etc...
    final_answer: str, output from the dataframe agent, which should contain a list of column names to be used in generating the view
    """
    prompt_RAG = (
        """
Context: {context}

The collection of Python scripts provided in the context, is designed to generate various types of data visualizations
using the mdvtools library. Each script focuses on a specific type of plot and follows a common structure that includes loading
data from a file, creating a plot using specific parameters, and serving the visualization through an MDV project.

All scripts in the context share a common workflow:

Setup: Define the project path, data path, and view name, the project path should always be: project_path = '"""
        + project.dir
        + """'
Plot function definition: Define the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) using a function in the same way as the context.
Project Creation: Initialize an MDVProject instance using the method: MDVProject(project_path, delete_existing=True).
Data Loading: Load data from the specified file into a pandas DataFrame using the load_data(data_path) function.
Data adding: Add the data source to the project using the method: project.add_datasource(datasource_name, data).
Plot Creation: Create the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) and define the plot paramaters in the same way as in the context.
Data Conversion: Convert the plot data to JSON format for integration with the MDV project using the convert_plot_to_json(plot) function.
Serving: Configure the project view, set it to editable, and serve the project using the .set_view(view_name, plot_view), .set_editable(True) and .serve() methods.

If no script is relevant, guided by the context generate a new script.

This string """
        + final_answer
        + """ specifies the names of the data fields that need to be plotted, for example in the params field. Get the structure of params definition from the context.

IMPORTANT: All string interpolations MUST use Python formatted string literals (f-strings).

Example:
# Correct:
message = f"gs|{{param}}(gs)|{{param_index}}"

# Common mistake to avoid:
WRONG: message = "gs|{{param}}(gs)|{{param_index}}"

The question for you to try to plot a graph based on the above instructions is given by this question: {question}""" + question + """

If the question asks for a gene-related query (e.g., gene expression value, most expressing gene, target expression, etc.), make sure that:
1. You load both datasources that are needed, e.g. cells and genes.
2. If gene expression is required, make sure the gene id, is given as a param in the formatted string literal format: f"gs|{{param}}(gs)|{{param_index}}", with param and param_index given by param = "param"
and param_index = data_frame_var.index.get_loc(param) respecitvely.

The data_path are given by this variable `""" + path_to_data + """`
The datasource_name is given by this variable `""" + datasouce_name + """`

"""
    )
    return prompt_RAG


def get_createproject_prompt_RAG_copy(project: MDVProject, path_to_data: str, datasource_name: str, final_answer: str, question: str) -> str:
    """
    Constructs a RAG prompt to guide LLM code generation for creating MDV plots.
    Handles both standard and gene-related queries.
    """
    prompt_RAG = (
        """
Context: {context}

The provided scripts demonstrate how to generate various data visualizations using the `mdvtools` library in Python.

Each script follows this standard workflow:

1. Setup:
    - Initialize an MDVProject instance using the method: MDVProject(project_path, delete_existing=True).
    - Use `scanpy.read_h5ad(data_path)` to load the AnnData object.

2. Data Loading:
    - Extract `adata.obs` into `data_frame_obs` (cell-level info).
    - Extract `adata.var` into `data_frame_var` (gene-level info).
    - Add a `name` column to `data_frame_var`: `adata.var_names.to_list()`

3. Datasource Registration:
    - Add data to the MDV project using:
        ```python
        project.add_datasource(datasource_name, data_frame_obs)
        project.add_datasource(datasource_name_2, data_frame_var)
        ```

4. Plot Construction:
    - Use a chart class (e.g., DotPlot, BoxPlot) and set `params = [...]` using selected fields.
        - The fields are given by """+final_answer+"""
    - Convert the chart to JSON using `convert_plot_to_json(plot)`
    - Set the view using `project.set_view(view_name, view_object)`

5. Parameter Handling:
    - The string """+final_answer+""" specifies the field names to use in the `params` list.
    - For parameters from `data_frame_obs` (cell-level), use them as-is.
    - For gene expression values from `data_frame_var`, use this syntax:
        ```python
        param = "GENE_NAME"
        param_index = data_frame_var['name'].tolist().index(param)
        f"gs|{{{{param}}}}(gs)|{{{{param_index}}}}"
        ```

6. Gene-Related Queries:
    If the question involves gene expression, expression comparison, or refers to gene names:
    - Load both `cells` and `genes` datasources.
    - Wrap gene names (from `data_frame_var`) using the syntax above.
    - Only wrap genes—do not apply `get_loc()` or `index` on `data_frame_obs` fields.

7. Your Task:
    - Interpret the user question and decide based on the question which graph needs to be plotted: """+question+"""
    - Use the fields """+final_answer+""" as params appropriately:
        - Wrap only gene names as shown.
        - Use others directly. Fields are case sensitive.
    - Use formatted f-strings for all dynamic strings.
    - Generate a valid Python script that creates and visualizes the appropriate chart using the MDVProject framework.
    - Update these variables with these values:
        - project_path = '"""+project.dir+"""'
        - data_path = '"""+path_to_data+"""'
        - view_name = a string (e.g., "default")
        - datasource_name = '"""+datasource_name+"""'

"""
    )
    return prompt_RAG