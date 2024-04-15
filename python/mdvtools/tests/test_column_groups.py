from mdvtools.mdvproject import MDVProject
import os
import pandas as pd

# make a project in local temp directory, which is ignored by git
path = os.path.join(os.path.dirname(__file__), "temp", "test_column_groups")
if not os.path.exists(path):
    os.makedirs(path)
p = MDVProject(path, delete_existing=True)

df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6], "c": [7, 8, 9]})

p.add_datasource("test", df)
