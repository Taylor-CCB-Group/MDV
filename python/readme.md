## Python tools for creating, manipulating and viewing MDV projects. 


### Installation

It is recommended, but not essential to create a virtual environment so there are no conflicts with modules in the global python.

To create and activate an environment in linux:-
```
python -m venv /path/to/myenv
source /path/to/myenv/bin/activate
````
In windows:-
```
python -m env c:\path\to\myenv
c:\path\to\myenv\Scripts\activate.bat
```

Installing mdv (with -e if for development) will include dependencies:
```
cd MDV/python
pip install -e .
```


### Quick Start

Then create an MDV project by giving the location of a folder
```python
from mdv.mdvproject import MDVProject
p = MDVProject("/path/to/hyp_example_data")
p.serve()
```
