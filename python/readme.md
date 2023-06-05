# Python tools for creating, manipulating and viewing MDV projects. 


## Installation

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

Next install the dependents 
```
pip install -r requirements.txt
```


## Quick Start

### Creating a Project

To make an empty project, just create an MDVProject object specifying a path. Data can then be added using the `add_datasource` method and specifying a pandas dataframe or the path to a text file.

```python
from mdv.mdvproject import MDVProject
p=MDVProject("/path/to/myproject")
df = pd.read_csv("cells.tsv",sep="\t")
p.add_datasource("cells",df)  #or p.add_datasource("cells","cells.tsv")
```

Initially an empty default view is created. The easiest way to add content is interactively through the browser. To do this, make the project editable and then serve it:-

```python
p.set_editable(True)
p.serve()
```

Then you can add charts and create views in the browser. Make sure that you save any changes you make.

### Serving the Project as a Static Page

To enable the project to be displayed as a static page without any backend architecture, the data needs to be converted to static binary files and the associated assets (JavaScript and images) added. To do this run the following:-
```python
p.convert_to_static_page(path/to/mywebproject)
```

You can then copy the mywebproject folder to the relevant location on your server, where content is served from. The project can be accessed with:-
```
http://mydomian.com/path/to/mywebproject
```
The index.html file contains a single div which houses the app. It can easily be altered to add custom headers/footers or other content e.g.
```html
    <div> My Header </div>
    <div id="holder"></div>
    <div> My Footer </div>
```

Because the JavaScript uses WebWorkers and SharedArrayBuffers, a [service worker](https://github.com/gzuidhof/coi-serviceworker) is added to load the necessary headers:- 
```
"Cross-Origin-Opener-Policy":"same-origin",
"Cross-Origin-Embedder-Policy":"require-corp"
```
If your server already adds these headers then the `convert_to_static_page` function can be invoked with `include_sab_headers=False`, which will exclude the service worker.


### Serving the Project in Development Mode
Create a static version of the project with the `debug=True` argument.
```
p.convert_to_static_page(/path/to/myproject, debug=True)
```
Then run webpack with the development config and point towards the folder you have just created.

```
webpack serve  --config dev.config.js --env folder=/path/to/myproject
```



