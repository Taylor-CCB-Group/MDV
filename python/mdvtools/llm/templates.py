
prompt_data = """
Based on the question asked and the dataframes provided, please perform the following steps:

1. Identify the data asked for in the question.
2. Based on step 1, find the relevant column names in the dataframes provided based on the information identified earlier in the question asked regarding data.
3. The relevant column names can be a combination from the two dataframes.
4. Provide the relevant column names from step 2 in a list.
5. If the question asks for a gene name, do provide a gene name in addition to the relevant columns names.
"""


packages_functions = """import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
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
from mdvtools.charts.violin_plot import ViolinPlot
from mdvtools.charts.multi_line_plot import MultiLinePlot
from mdvtools.charts.table_plot import TablePlot
from mdvtools.charts.wordcloud_plot import WordcloudPlot
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

# this is Maria's un-sullied prompt - apart from the project name
def get_createproject_prompt_RAG(project_id: str, path_to_data: str, final_answer: str):
    prompt_RAG = (
""" 
Context: {context}

The collection of Python scripts provided in the context, is designed to generate various types of data visualizations 
using the mdvtools library. Each script focuses on a specific type of plot and follows a common structure that includes loading 
data from a file, creating a plot using specific parameters, and serving the visualization through an MDV project. 

All scripts in the context share a common workflow:

Setup: Define the project path, data path, and view name, the project path should always be: project_path = os.path.expanduser('~/mdv/project')
Plot function definition: Define the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) using a function in the same way as the context.
Project Creation: Initialize an MDVProject instance using the method: MDVProject(project_path, delete_existing=True).
Data Loading: Load data from the specified file into a pandas DataFrame using the load_data(path) function.
Data adding: Add the data source to the project using the method: project.add_datasource(data_path, data).
Plot Creation: Create the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) and define the plot paramaters in the same way as in the context.
Data Conversion: Convert the plot data to JSON format for integration with the MDV project using the convert_plot_to_json(plot) function.
Serving: Configure the project view, set it to editable, and serve the project using the .set_view(view_name, plot_view), .set_editable(True) and .serve() methods.

You are a top-class Python developer. Based on the question: {question}, decide which script from the context {context} is more relevant to the question: {question} and update the script to address the question.
If no script is relevant, guided by the context generate a new script. 

This list """ + final_answer + """ specifies the names of the data fields that need to be plotted, for example in the params field. Get the structure of params definition from the context.
DO NOT forget to use the f-string, or formatted string literal, python structure in the parameters, params or param.

If the prompt asks for a gene, make sure you load this datasource and that you create a link between the two datasets.

The path to data are given by this variable """ + path_to_data + """
"""
    )
    return prompt_RAG

# The data should be loaded in the same way as in the examples. 

# different prompt for experimenting with update project
def get_updateproject_prompt_RAG(project_name: str, path_to_data: str, final_answer: str):
    prompt_RAG = (
""" You can only generate python code based on the provided context. If a response cannot be formed strictly using the context, politely say you need more information to generate the plot.

Context: {context}]

The collection of Python scripts provided in the context, is designed to generate various types of data visualizations 
using the mdvtools library. Each script focuses on a specific type of plot and follows a common structure that includes loading 
an existing MDV project, creating a plot using specific parameters, and adding the visualization to the project. 

All scripts in the context share a common workflow:

Setup: Define the project path, data path, and view name, the project path should always be: project_path = os.path.expanduser('~/mdv/""" + project_name + """')
Plot function definition: Define the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) using a function in the same way as the context.
Project Creation: Initialize an MDVProject instance using the method: MDVProject(project_path).
Data Loading: Load data from the specified CSV file into a pandas DataFrame using the load_data(path) function.
Data adding: Add the data source to the project using the method: project.add_datasource(data_path, data).
Plot Creation: Create the respective plot (dot plot, heatmap, histogram, box plot, scatter plot, 3D scatter plot, pie/ring chart, stacked row plot) and define the plot paramaters in the same way as in the context.
Data Conversion: Convert the plot data to JSON format for integration with the MDV project using the convert_plot_to_json(plot) function.

You are a top-class Python developer. Generate a Python script following the workflow detailed above and use exactly the same lines of code as the scripts in the context. 
The data should be loaded from a csv provided, the path to the csv is given by: """ + path_to_data + """
This list """ + final_answer + """ specifies the names of the data fields that need to be plotted, for example in the params field. Get the structure of params definition from the context. 
The question: {question} specifies the plot required. 
"""
    )
    return prompt_RAG