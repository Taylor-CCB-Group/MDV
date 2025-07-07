from mdvtools.mdvproject import MDVProject

def create_project_markdown(project: MDVProject) -> str:
    markdown = f"\n\n**Available Columns:**\n\n"
    for name in project.get_datasource_names():
        ds = project.get_datasource_metadata(name)
        markdown += f"**{name}:** ({ds['size']} rows)\n\n"
        columns = [(c['name'], c['datatype']) for c in ds['columns']]
        markdown += create_column_markdown(columns)

    return markdown

def create_column_markdown(cols: list[tuple[str, str]]) -> str:
    if cols:
        markdown = "\n| Column Name | Data Type |\n|---|---|\n"
        for col in cols:
            markdown += f"| {col[0]} | {col[1]} |\n"
    else:
        markdown = "No columns found."
    return f"{markdown}\n\n"