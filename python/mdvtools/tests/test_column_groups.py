from mdvtools.mdvproject import MDVProject
import os
import pandas as pd


def test_add_datasource_with_multiple_columns():
    # Make a project in a local temp directory that is ignored by git.
    path = os.path.join(os.path.dirname(__file__), "temp", "test_column_groups")
    os.makedirs(path, exist_ok=True)

    project = MDVProject(path, delete_existing=True)
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})
    project.add_datasource("test", df)

    datasource_columns = project.datasources[0]["columns"]
    assert [column["name"] for column in datasource_columns] == ["a", "b", "c"]
