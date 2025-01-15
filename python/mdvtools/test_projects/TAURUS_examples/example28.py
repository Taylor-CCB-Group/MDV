import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.pie_chart import PieChart
import json 




def load_data(path):
    #Load data from the specified CSV file.
    return pd.read_csv(path, low_memory=False)

def convert_plot_to_json(plot):
    #Convert plot data to JSON format.
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))
    

def create_pie_chart(title, param, size, position):
    plot = PieChart(
        title=title,
        param=param,
        size=size,
        position=position
    )
    return plot

def convert_plot_to_json(plot):
    return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\", ""))

def main():
    project_path = os.path.expanduser('~/mdv/project')
    view_name = "TAURUS Cell Type Composition"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    cells_df.name = 'cells'
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource('cells', cells_df)
    
    # PieChart parameters
    pie_title = "Cell Type Composition in TAURUS Dataset"
    pie_params = "final_analysis"  # Assuming 'final_analysis' is the column for cell types
    pie_size = [400, 400]
    pie_position = [10, 10]
    
    # Create pie chart
    pie_chart = create_pie_chart(pie_title, pie_params, pie_size, pie_position)
    
    # Convert plot to JSON and set view
    pie_chart_json = convert_plot_to_json(pie_chart)
    pie_chart_view = {'initialCharts': {'cells': [pie_chart_json]}}
    
    project.set_view(view_name, pie_chart_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()