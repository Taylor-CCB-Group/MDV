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

Next install the dependencies
```
pip install -r requirements.txt
```


### Quick Start - to serve an HDF5 based project from a local python server

Then create an MDV project by giving the location of a folder
```python
from mdv.mdvproject import MDVProject
p = MDVProject("/path/to/hyp_example_data")
p.serve()
```

### Quick Start - to convert CSV data for a simple file server

Given a CSV file associated folder of images, the `csv_to_static` script will create a folder (containing `datasources.json`, `state.json`, `views.json`, `{mydata.csv}.b`, `{mydata.csv}.json`, `images/`) that can be served by a simple file server.

The script takes arguments for the CSV file and the output folder. If either of these is not present, or if the input file is not found, it will prompt for these interactively.

```
(myenv) $ python csv_to_static.py -i /path/to/mydata.csv -o /path/to/myoutputfolder
```

#### Example Python code to serve the static data

The following code will serve the static data from a local python server. This is useful for testing, but not recommended for production use. Running this code will serve the data on http://localhost:8000, such that it can then be viewed at https://mdv-dev.netlify.app/src/static.html?dir=http://localhost:8000

```python
import http.server
import socketserver


class CustomHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type, responsetype')
        self.send_header('Cross-Origin-Opener-Policy', 'cross-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
        self.send_header('Cross-Origin-Resource-Policy', 'cross-origin')
        super().end_headers()

    def do_OPTIONS(self):
        self.send_response(200)
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type, responsetype')
        self.send_header('Cross-Origin-Opener-Policy', 'cross-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
        self.send_header('Cross-Origin-Resource-Policy', 'cross-origin')
        super().end_headers()

with socketserver.TCPServer(("", 8000), CustomHandler) as httpd:
    print("serving at port", 8000)
    httpd.serve_forever()
```