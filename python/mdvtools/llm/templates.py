from mdvtools.mdvproject import MDVProject
from typing import Any
prompt_data = """
Your task is to:  
1. Identify the type of data the user needs (e.g., categorical, numerical, etc.) by inspecting the DataFrames provided.
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
    - Sankey plot: Requires two categorical columns.  
      - If only one categorical variable is available, return it twice.  
    - Scatter plot (2D): Requires two numerical columns and one any column for color.
    - Scatter plot (3D): Requires three numerical columns and one any column for color.  
    - Selection dialog plot: Requires any column.  
    - Stacked row chart: Requires two categorical columns.  
      - If only one categorical variable is available, return it twice.  
    - Table Plot: Requires any column(s).  
    - Text box: Requires no columns, just text.  
    - Violin plot: Requires one categorical column and one numerical column.  
    - Wordcloud: Requires one categorical column.  
8. Important: Clearly separate the selected columns with quotes and commas.
9. The column names are case sensitive therefore return them as they are defined in the dataframe.
10. Output format:
   - First line: Quoted, comma-separated list of column names.
   - Second line: Quoted, comma-separated list of suitable chart types for the selected columns.
11. NEVER explain your reasoning.
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



def get_createproject_prompt_RAG(project: MDVProject, path_to_data: str, datasource_name: str, final_answer: str, question: str) -> str:
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
    - Use a chart class (e.g., DotPlot, BoxPlot, SelectionDialogPlot) and set `params = [...]` using selected fields.
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

7. Queries requiring subsetting of the dataset:
    If to answer the question requires a subset of the data or filtering the data, make sure to:
    - Add a selection dialog plot with all the parameters that were passed on as params.

8. Your Task:
    - Interpret the user question and decide based on the question which graph needs to be plotted: """+question+"""
    - Use the fields """+final_answer+""" as params appropriately:
        - Wrap only gene names as shown.
        - Use others directly. Fields are case sensitive.
    - Use formatted f-strings for all dynamic strings.
    - Generate a valid Python script that creates and visualizes the appropriate chart using the MDVProject framework.
    - Update these variables with these values:
        - project_path = '"""+project.dir+"""'
        - data_path = '"""+path_to_data+"""'
        - view_name = a string describing what is being visualized.
        - datasource_name = '"""+datasource_name+"""'
    - The type of fields given in params should follow the following visualisation guidelines for each type of chart:
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
        - Sankey plot: Requires two categorical columns.  
            - If only one categorical variable is available, return it twice.  
        - Scatter plot (2D): Requires two numerical columns and one any column for color.
        - Scatter plot (3D): Requires three numerical columns and one any column for color.  
        - Selection dialog plot: Requires any column.  
        - Stacked row chart: Requires two categorical columns.  
            - If only one categorical variable is available, return it twice.  
        - Table Plot: Requires any column(s).  
        - Text box: Requires no columns, just text.  
        - Violin plot: Requires one categorical column and one numerical column.  
        - Wordcloud: Requires one categorical column.
"""
    )
    return prompt_RAG
