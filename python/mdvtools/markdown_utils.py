from mdvtools.mdvproject import MDVProject
from typing import Optional
import json

def create_error_markdown(message: str, traceback: Optional[str] = None, extra_metadata: Optional[dict] = None) -> str:
    """
    Create a markdown representation of an error.
    Args:
        message: The error message
        traceback: The traceback of the error
        extra_metadata: Extra metadata about the error
    Returns:
        A string containing the markdown representation of the error
    """
    markdown = f"**Error:** {escape_markdown(message)}\n\n"
    if traceback:
        # Use HTML details tag for collapsable section
        markdown += f"<details><summary>Traceback</summary>\n\n```\n{traceback}\n```\n\n</details>\n\n"
    if extra_metadata:
        markdown += f"<details><summary>Extra Metadata</summary>\n\n```json\n{json.dumps(extra_metadata, indent=2)}\n```\n\n</details>\n\n"
    return markdown

def create_project_markdown(project: MDVProject, wrap_in_details: bool = True) -> str:
    """
    Create a markdown representation of the project.
    Args:
        project: The MDVProject object
        wrap_in_details: Whether to wrap the markdown in a details tag
    Returns:
        A string containing the markdown representation of the project
    """
    markdown = ""
    if wrap_in_details:
        markdown += "\n\n<details>\n\n"
        markdown += f"<summary>Project Datasources</summary>\n\n"
    for name in project.get_datasource_names():
        ds = project.get_datasource_metadata(name)
        markdown += f"## **{name}:** ({ds['size']} rows)\n\n"
        # todo - add a summary of the data here
        markdown += create_column_markdown(ds["columns"])
    if wrap_in_details:
        markdown += "\n\n</details>\n\n"
    return markdown

def escape_markdown(text: str) -> str:
    """
    Escape special markdown characters.
    In future we may want more advanced markdown functionality which may warrant
    adding a markdown library.
    """
    special_chars = "\\`*_{}[]()#+-.!<>|"
    for char in special_chars:
        text = text.replace(char, "\\" + char)
    return text

def create_column_markdown(cols: list[dict]) -> str:
    """
    Create a markdown representation of the columns in a datasource.
    Args:
        cols: A list of dictionaries, each containing the column name and data type
    Returns:
        A string containing the markdown representation of the columns
    """
    numeric_cols = []
    categorical_cols = []
    other_cols = []

    for c in cols:
        dt = c.get('datatype')
        if dt in ['integer', 'double', 'int32']:
            numeric_cols.append(c)
        elif dt in ['text', 'text16', 'multitext', 'unique']:
            categorical_cols.append(c)
        else:
            other_cols.append(c)

    markdown = ""
    
    if categorical_cols:
        markdown += "\n### Categorical Columns\n"
        markdown += "| Column Name | Data Type | Values (sample) |\n|---|---|---|\n"
        for col in categorical_cols:
            col_name = escape_markdown(col.get('name', ''))
            datatype = col.get('datatype', '')
            values_str = ""
            if col.get('values'):
                sample_values = col['values'][:5]
                escaped_values = [escape_markdown(str(v)) for v in sample_values]
                values_str = ", ".join(escaped_values)
                if len(col['values']) > 5:
                    values_str += ", ..."
            markdown += f"| {col_name} | {datatype} | {values_str} |\n"

    if numeric_cols:
        markdown += "\n### Numeric Columns\n"
        markdown += "| Column Name | Data Type | Min / Max | Quantiles (0.05) |\n|---|---|---|---|\n"
        for col in numeric_cols:
            col_name = escape_markdown(col.get('name', ''))
            datatype = col.get('datatype', '')
            min_max = col.get('minMax')
            min_max_str = f"{min_max[0]} / {min_max[1]}" if min_max and len(min_max) == 2 else ""
            
            quantiles = col.get('quantiles', {})
            q05 = quantiles.get('0.05')
            quantiles_str = f"[{q05[0]}, {q05[1]}]" if q05 else ""

            markdown += f"| {col_name} | {datatype} | {min_max_str} | {quantiles_str} |\n"

    if other_cols:
        markdown += "\n### Other Columns\n"
        markdown += "| Column Name | Data Type |\n|---|---|\n"
        for col in other_cols:
            col_name = escape_markdown(col.get('name', ''))
            datatype = col.get('datatype', 'N/A')
            markdown += f"| {col_name} | {datatype} |\n"
    
    if not (numeric_cols or categorical_cols or other_cols):
        markdown = "No columns found."
    
    return f"{markdown}\n\n"

## ---------------------------------
## ---------------------------------
## chat/LLM related stuff... there could be an argument for having this in `mdvtools.llm`
## but it doesn't have any external module dependencies etc.


chart_types_md = """
- Abundance Box Plot
- Box Plot
- Density Scatter Plot
- Dot Plot
- Heat Map
- Histogram Plot
- Multi Line Chart
- Pie Chart
- Row Chart
- Row Summary Box
- Sankey Diagram
- Stacked Row Chart (Categorical Heatmap)
- Table
- Text Box
- 2D Scatter Plot
- Violin Plot
- Word Cloud
"""

example_intents_md = """
- "distribution", "spread" = Histogram Plot, Box Plot, Violin Plot
- "relationship", "correlation" = Scatter Plot, Density Scatter, Heat Map
- "comparison", "difference", "change" = Box Plot, Violin Plot, Multi Line Chart, Dot Plot
- "composition", "proportion", "breakdown" = Pie Chart, Stacked Row Chart, Row Chart
- "over time", "trend", "temporal" = Line Chart, Multi Line Chart
- "expression", "gene", "marker" = Dot Plot, Heat Map, Box Plot
- "spatial", "location", "embedding" = 2D Scatter Plot, Density Scatter Plot
- "flow", "transition" = Sankey Diagram
- "metadata", "category", "annotation" = Table, Row Summary Box, Row Chart
- "filter", "subset", "select" = Selection Dialog Plot
"""

def create_suggested_questions_prompt(project: MDVProject) -> str:
    return f"""
    You are an expert in the data in this project.
    You are given a project with the following datasources:
    {create_project_markdown(project)}

    Given this data, generate a list of 5 biologically relevant questions we could visualise.
    When phrasing the questions, only refer to column names found in that project data, and only suggest visualisations using these charts:

    {chart_types_md}

    ## Example intents:

    {example_intents_md}

    If a query needs multiple charts to be plotted, show them as individual queries. 
    For example: What is the distribution of age_or_mean_of_age_range or BMI among different disease categories? 
    Chart: Histogram Plot, Box Plot, Violin Plot should be What is the distribution of age_or_mean_of_age_range 
    among different disease categories? Violin Plot to keep the visualisation simple.

    Each question should be a single, complete sentence.
    """