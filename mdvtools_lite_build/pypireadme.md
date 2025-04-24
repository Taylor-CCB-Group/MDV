# MDV Tools Lite

This is a lightweight version of MDV, designed to simplify project creation and deployment. It provides:

- Tools for creating and manipulating MDV projects
- A **lightweight server** for viewing and modifying projects 
- A method to  **create static websites**, making it ideal for projects that require minimal server-side dependencies.


To install from  TestPyPi:-
```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ mdvtools
```

To test, run this script (change /path/to/project to something appropriate):

```python
import scanpy
from mdvtools.conversions import convert_scanpy_to_mdv
folder = "/path/to/mytestproject"
data = scanpy.datasets.pbmc3k_processed()
proj= convert_scanpy_to_mdv(folder,data)
proj.serve(port=5057)
```

The MDV project will be in the folder you specified and you can view it by pointing your browser to localhost:5057


