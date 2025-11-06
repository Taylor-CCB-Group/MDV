import os
import pandas as pd
import scanpy as sc
from mdvtools.mdvproject import MDVProject
from mdvtools.charts.pie_chart import PieChart
import json 

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
    view_name = "default"
    
    # Load data
    data_path = "file_path"
    adata = sc.read_h5ad(data_path)
    cells_df = pd.DataFrame(adata.obs)
    datasource_name = "datasource_name"
    
    # Filter data for sample A (assuming sample A is 'CID003352-2')
    sample_a_df = cells_df[cells_df['sample_id'] == 'CID003352-2']
    assert isinstance(sample_a_df, pd.DataFrame)
    
    # Create project
    project = MDVProject(project_path, delete_existing=True)
    
    # Add datasource
    project.add_datasource(datasource_name, sample_a_df)
    
    # PieChart parameters
    pie_title = "Cell Type Proportions in Sample A"
    pie_params = "final_analysis"  # Assuming 'final_analysis' is the column for cell types
    pie_size = [400, 400]
    pie_position = [10, 10]
    
    # Create pie chart
    pie_chart = create_pie_chart(pie_title, pie_params, pie_size, pie_position)
    
    # Convert plot to JSON and set view
    pie_chart_json = convert_plot_to_json(pie_chart)
    pie_chart_view = {'initialCharts': {datasource_name: [pie_chart_json]}}
    
    project.set_view(view_name, pie_chart_view)
    project.set_editable(True)
    project.serve()

if __name__ == "__main__":
    main()