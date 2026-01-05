from mdvtools.mdvproject import MDVProject
from typing import Any
from mdvtools.llm.markdown_utils import create_project_markdown, create_column_markdown
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
           - If it does not exist, assume the user provided a gene name and that df2["name"] may contain gene IDs instead.
                - Attempt to match the user-provided gene name to the corresponding gene ID using any available mapping logic (e.g., a lookup function or mapping dictionary).
                - If a corresponding gene ID is found in df2["name"], return that value.
           - If it exists, return it.
           - If no match is found, ignore the requested gene and instead select one or more gene names from df2["name"].
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
    - Box plot: Requires only one categorical column and one numerical column.  
    - Density Scatter plot: Requires two numerical columns and one categorical column.  
    - Dot plot: Requires only one categorical column and any number of numerical columns.  
    - Heatmap: Requires only one categorical column and any number of numerical columns.  
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
    - Violin plot: Requires only one categorical column and one numerical column.  
    - Wordcloud: Requires one categorical column.  
8. Important: Clearly separate the selected columns with quotes and commas.
9. The column names are case sensitive therefore return them as they are defined in the dataframe.
10. Output format:
   - First line: The word "fields" following by quoted, comma-separated list of column names.
   - Second line: The word "charts" following by quoted, comma-separated list of suitable chart types for the selected columns.
11. NEVER explain your reasoning.
"""

packages_functions = """import os
import pandas as pd
try:
    import scanpy as sc
    HAS_SCANPY = True
except Exception:
    sc = None
    HAS_SCANPY = False
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

# Ensure all charts are available in the local namespace even if ruff tries to be clever
_ = [MDVProject, DensityScatterPlot, HeatmapPlot, HistogramPlot, DotPlot, BoxPlot, ScatterPlot3D, RowChart, ScatterPlot, AbundanceBoxPlot, StackedRowChart, RingChart, PieChart, ViolinPlot, MultiLinePlot, TablePlot, WordcloudPlot, TextBox, RowSummaryBox, SelectionDialogPlot, SankeyPlot]
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
    # Build markdown context for the selected datasource; fall back to whole project if needed
    try:
        ds_meta = project.get_datasource_metadata(datasource_name)
        context_md = f"## **{datasource_name}:** ({ds_meta['size']} rows)\n\n" + create_column_markdown(ds_meta["columns"])
    except Exception:
        context_md = create_project_markdown(project)
    prompt_RAG = (
    """
    Project Data Context:

    """ + context_md + """

    Context: {context}

    The provided scripts demonstrate how to generate various data visualizations using the `mdvtools` library in Python.

    Each script follows this workflow:

    1. Setup:
        - Initialize an MDVProject instance using the method: MDVProject(project_path, delete_existing=True) when creating a new project.
        - When modifying an existing project (the default in this chat context), do NOT recreate or delete the project.

    2. Data Access (scanpy optional):
        - If HAS_SCANPY is True and `data_path` points to a valid .h5ad file, you MAY load AnnData:
            adata = sc.read_h5ad(data_path)
            data_frame_obs = adata.obs
            data_frame_var = adata.var.assign(name=adata.var_names.to_list())
        - OTHERWISE (no scanpy or no .h5ad):
            - DO NOT use scanpy.
            - Use the existing MDV datasources already registered in the project:
                data_frame_obs = project.get_datasource_as_dataframe('cells') if available
                data_frame_var = project.get_datasource_as_dataframe('genes') if available
            - If a 'genes' datasource is not available, omit gene-specific logic and charts that require gene wrapping.

    3. Datasource Registration:
        - When modifying an existing project, DO NOT call project.add_datasource(...).
        - Only add datasources when explicitly creating a new project from raw data (not the default).

    4. Plot Construction:
        - Use a chart class (e.g., DotPlot, BoxPlot, SelectionDialogPlot) and set `params = [...]` using selected fields.
            - The fields and chart type are given by """+final_answer+"""
        - Convert the chart to JSON using `convert_plot_to_json(plot)`
        - Set the view using `project.set_view(view_name, view_object)`

    5. Parameter Handling:
        - The string """+final_answer+""" specifies the field names to use in the `params` list and the chart type to use.
        - For parameters from cell-level data (from data_frame_obs or existing 'cells' datasource), use them as-is.
        - For gene expression (if a 'genes' datasource or data_frame_var is available), use this syntax to refer to a gene:
            ```python
            param = "GENE_NAME"
            param_index = data_frame_var['name'].tolist().index(param)
            f"gs|{{param}}(gs)|{{param_index}}"
            ```

    6. Gene-Related Queries:
        - Only perform gene wrapping if you have a usable gene names table (data_frame_var or 'genes' datasource with a 'name' column).
        - If not available, proceed without gene wrapping and prefer non-gene charts.

    7. Queries requiring subsetting of the dataset:
        If to answer the question requires a subset of the data or filtering the data, make sure to:
        - Add a selection dialog plot with all the parameters that were passed on as params. Make sure it has a title.

    8. Always add a selection dialog plot along the other charts. It must have all the parameters that were passed on as params. It must have a title.

    9. Your Task:
        - Interpret the user question and decide based on the question which graph needs to be plotted: """+question+final_answer+"""
        - Use the fields in the """+final_answer+""" as params appropriately:
            - Wrap only gene names as shown.
            - Use others directly. Fields are case sensitive.
        - Use formatted f-strings for all dynamic strings.
        - Generate a valid Python script that creates and visualizes the appropriate chart using the MDVProject framework.
        - Update these variables with these values:
            - project_path = '"""+project.dir+"""'
            - data_path = '"""+path_to_data+"""'  # may be empty; if empty or HAS_SCANPY is False, DO NOT use scanpy.
            - view_name = a string, in double quotes, describing what is being visualized.
            - datasource_name = '"""+datasource_name+"""'
        - The possible charts are given by """+final_answer+""" and should follow the following visualisation guidelines for each type of chart:
            - Abundance Box plot: Requires three categorical columns.  
                - If only one categorical variable is available, use it three times.  
                - If two are available, use one of them twice.  
            - Box plot: Requires only one categorical column and one numerical column.  
            - Density Scatter plot: Requires two numerical columns and one categorical column.  
            - Dot plot: Requires only one categorical column and any number of numerical columns.  
            - Heatmap: Requires only one categorical column and any number of numerical columns.  
            - Histogram: Requires one numerical column.  
            - Multiline chart: Requires one numerical column and one categorical column.  
            - Pie Chart: Requires one categorical column.  
            - Row Chart: Requires one categorical column.  
            - Row summary box: Requires any column(s).  
            - Sankey plot: Requires two categorical columns.  
                - If only one categorical variable is available, use it twice.  
            - Scatter plot (2D): Requires two numerical columns and one any column for color.
            - Scatter plot (3D): Requires three numerical columns and one any column for color.  
            - Selection dialog plot: Requires any column.  
            - Stacked row chart: Requires two categorical columns.  
                - If only one categorical variable is available, use it twice.  
            - Table Plot: Requires any column(s).  
            - Text box plot: Requires no columns, just text.  
            - Violin plot: Requires only one categorical column and one numerical column.  
            - Wordcloud: Requires one categorical column.
    Output format: Return the python code that is to be run to generate the charts.

    After generating the code, include a detailed explanation in your response (renderable by markdown renderer) that covers:
    1. Why this chart is the best way to answer the question
    2. What biological insights can be gained from this visualization
    3. What subsequent analysis tasks could be performed based on these results

    The explanation will be displayed in the chat window automatically.



"""
    )
    return prompt_RAG
