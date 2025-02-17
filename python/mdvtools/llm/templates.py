from mdvtools.mdvproject import MDVProject

prompt_data = """
Your task is to:  
1. Identify the type of data the user needs (e.g., categorical, numerical, etc.).  
2. Select the most relevant column names from the first DataFrame provided unless handling a gene-related query.  
3. If the query is gene-related (e.g., most expressing gene, target expression, etc.), retrieve gene names 
from the second DataFrame while selecting the remaining columns from the first DataFrame.  
4. Ensure that the selected columns match the visualization requirements:  
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
    - Row summary box: Requires any column.  
    - Sankey diagram: Requires two categorical columns.  
      - If only one categorical variable is available, return it twice.  
    - Scatter plot (2D): Requires two numerical columns.  
    - Scatter plot (3D): Requires three numerical columns and one categorical column for color.  
    - Selection dialog plot: Requires any column.  
    - Stacked row chart: Requires two categorical columns.  
      - If only one categorical variable is available, return it twice.  
    - Table: Requires any column.  
    - Text box: Requires no columns, just text.  
    - Violin plot: Requires one categorical column and one numerical column.  
    - Wordcloud: Requires one categorical column.  
5. Return the column names in a string format, e.g., "col1", "col2".  
6. Do not provide additional explanations—only return the string.
"""


packages_functions = """import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
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

def load_data(path):
    #Load data from the specified CSV file.
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    #Convert plot data to JSON format.
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))

"""

def get_createproject_prompt_RAG(project: MDVProject, path_to_data: str, datasouce_name: any, final_answer: str) -> str:
    """
    Returns the prompt for the create project RAG.
    
    Args:
    project: MDVProject we pass in the project object to get the project directory - previously we were passing the project id
      and making (wrong) assumptions about the project directory based on that - while also being less.
    path_to_data: str, path to an anndata object, subject to change
    datasouce_name: any, we need to better clarify how we reason about multi-ds queries etc...
    final_answer: str, output from the dataframe agent, which should contain a list of column names to be used in generating the view
        - note: what does this mean in cases where we are trying to make multiple charts, belonging to different datasources, etc?
    """
    # assert isinstance(project, MDVProject) # this has problems, with autoreload? was erroneously failing...
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

You are a top-class Python developer. Based on the question: {question}, decide which script from the context {context} is more relevant to the question: {question} and update the script to address the question.
If no script is relevant, guided by the context generate a new script.

This list """
        + final_answer
        + """ specifies the names of the data fields that need to be plotted, for example in the params field. Get the structure of params definition from the context.
DO NOT forget to use the f-string, or formatted string literal, python structure in the parameters, params or param.

Generate Python code that correctly uses formatted string literals.
Ensure all strings requiring variable substitution use the 'f' prefix.

Example:
# Correct:
message = f"Hello, {question}!"

Incorrect:
message = "Hello, {question}!"

If the prompt asks for a gene, make sure you load this datasource and that you create a link between the two datasets.

The data_path are given by this variable `""" + path_to_data + """`
The datasource_name is given by this variable `""" + datasouce_name + """`
"""
    )
    return prompt_RAG

# The data should be loaded in the same way as in the examples.

# different prompt for experimenting with update project - unused as of now
def get_updateproject_prompt_RAG(project_name: str, ds_name: str, final_answer: str):
    prompt_RAG = (
        """
Context: {context}]

The collection of Python scripts provided in the context, is designed to generate various types of data visualizations
using the mdvtools library. Each script focuses on a specific type of plot and follows a common structure that includes loading
an existing MDV project, creating a plot using specific parameters, and adding the visualization to the project.

All scripts in the context share a common workflow:

Setup: Define the project path, which should always be: project_path = os.path.expanduser('~/mdv/"""
        + project_name
        + """')
Plot function definition: Define the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) using a function in the same way as the context.
Project Creation: Initialize an MDVProject instance using the method: MDVProject(project_path).
Data Loading: Data can be retrieved from an existing datasource in the project as a pandas DataFrame using the method: project.get_datasource_as_dataframe(ds_name).
Plot Creation: Create the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) and define the plot paramaters in the same way as in the context.
Data Conversion: Convert the plot data to JSON format for integration with the MDV project using the convert_plot_to_json(plot) function.

You are a top-class Python developer. Based on the question: {question}, decide which script from the context {context} is more relevant to the question: {question} and update the script to address the question.
If no script is relevant, guided by the context generate a new script.

This list """
        + final_answer
        + """ specifies the names of the data fields that need to be plotted, for example in the params field. Get the structure of params definition from the context.
DO NOT forget to use the f-string, or formatted string literal, python structure in the parameters, params or param.
"""
    )
    return prompt_RAG
