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
        markdown += f"**{name}:** ({ds['size']} rows)\n\n"
        columns = [(c['name'], c['datatype']) for c in ds['columns']]
        markdown += create_column_markdown(columns)

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

def create_column_markdown(cols: list[tuple[str, str]]) -> str:
    if cols:
        markdown = "\n| Column Name | Data Type |\n|---|---|\n"
        for col in cols:
            # todo - consider more comprehensive escaping and validation
            # particularly if we use this for e.g. LLM queries.
            col_name = escape_markdown(col[0])
            markdown += f"| {col_name} | {col[1]} |\n"
    else:
        markdown = "No columns found."
    return f"{markdown}\n\n"