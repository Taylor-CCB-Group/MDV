from mdvtools.mdvproject import MDVProject
from typing import Optional
import json

def create_error_markdown(message: str, traceback: Optional[str] = None, extra_metadata: Optional[dict] = None) -> str:
    markdown = f"**Error:** {escape_markdown(message)}\n\n"
    if traceback:
        # Use HTML details tag for collapsable section
        markdown += f"<details><summary>Traceback</summary>\n\n```\n{traceback}\n```\n\n</details>\n\n"
    if extra_metadata:
        markdown += f"<details><summary>Extra Metadata</summary>\n\n```json\n{json.dumps(extra_metadata, indent=2)}\n```\n\n</details>\n\n"
    return markdown

def create_project_markdown(project: MDVProject) -> str:
    markdown = "\n\n<details>\n\n"
    markdown += f"<summary>Project Datasources</summary>\n\n"
    for name in project.get_datasource_names():
        ds = project.get_datasource_metadata(name)
        markdown += f"## **{name}:** ({ds['size']} rows)\n\n"
        # todo - add a summary of the data here
        markdown += create_column_markdown(ds["columns"])

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