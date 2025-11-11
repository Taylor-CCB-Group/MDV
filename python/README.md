# Python tools for creating, manipulating and viewing MDV projects. 


## Installation

### Make virtual environment

It is recommended, but not essential to create a virtual environment so there are no conflicts with modules in the global python.

To create and activate an environment in a Unix-like system:-
```
python -m venv /path/to/myenv
source /path/to/myenv/bin/activate
````
In windows:-
```
python -m venv c:\path\to\myenv
c:\path\to\myenv\Scripts\activate.bat
```

### Install poetry

Install poetry if it is not already installed. This can be done with the official installer:

```
curl -sSL https://install.python-poetry.org | python3 -
```

Or with pipx:
```
pip install pipx
pipx install poetry
```

See the [poetry installation instructions](https://python-poetry.org/docs/#installing-with-pipx) for more details and troubleshooting.

### Install MDV

To install MDV, run:

```
cd MDV/python
poetry install --with dev
```

## Quick Start

### Creating a Project

To make an empty project, just create an MDVProject object specifying a path. Data can then be added using the `add_datasource` method and specifying a pandas dataframe or the path to a text file.

### Quick Start - to serve an HDF5 based project from a local python server

```python
from mdvtools.mdvproject import MDVProject

p = MDVProject("/path/to/myproject")
df = pd.read_csv("cells.tsv", sep="\t")
p.add_datasource("cells", df)  # or p.add_datasource("cells", "cells.tsv")

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

You can then copy the mywebproject folder to the relevant location on your server, where content is served from. It needs to be served over https because a service worker is used to add headers. It will also need to be able to use range headers for loading portions of binary data. The project can be accessed with:-
```
https://mydomian.com/path/to/mywebproject
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
If your server already adds these headers then the `convert_to_static_page` function can be invoked with `include_sab_headers=False`, which will exclude the service worker. In this case the project does not have to be served over https

### To convert CSV data for a simple file server

**This is somewhat deprecated in favour of the `convert_to_static_page` method and will likely be removed.**

Given a CSV file associated folder of images, the `csv_to_static` script will create a folder (containing `datasources.json`, `state.json`, `views.json`, `{mydata.csv}.b`, `{mydata.csv}.json`, `images/`) that can be served by a simple file server. 

```

The script takes arguments for the CSV file and the output folder. If either of these is not present, or if the input file is not found, it will prompt for these interactively.

```
(myenv) $ python csv_to_static.py -i /path/to/mydata.csv -o /path/to/myoutputfolder
```

#### Example Python code to serve the static data

The following code will serve the static data from a local python server. This is useful for testing, but not recommended for production use. Running this code will serve the data on http://localhost:8000, such that it can then be viewed at https://mdv-dev.netlify.app/?dir=http://localhost:8000

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

## Configuration

### User Configuration File

MDV supports user-provided configuration for extensions via a JSON file specified through an environment variable:

- **MDV_USER_CONFIG_PATH**: Set this environment variable to the path of your user configuration file to configure extensions. Available extensions include:
  - `chat`: Enables the chat/LLM extension
  - `project_manager`: Enables the project management extension
  
  Example user_config.json:
  ```json
  {
    "extensions": ["chat", "project_manager"]
  }
  ```
  
  Example usage:
  ```bash
  export MDV_USER_CONFIG_PATH="/path/to/your/user_config.json"
  ```
  
  **Important for Docker deployments**: You must mount the config file into the container so it can access it. For example, if your config file is at `./config/user_config.json` on the host, you would:
  
  1. Set the environment variable: `MDV_USER_CONFIG_PATH=/config/user_config.json`
  2. Mount the file in docker-compose: `- ./config/user_config.json:/config/user_config.json`
  
  **Extension Configuration Behavior**:
  
  - If `MDV_USER_CONFIG_PATH` is **set** and the file exists and is valid: Extensions are loaded from the user config file.
  - If `MDV_USER_CONFIG_PATH` is **set** but the file is missing or cannot be parsed: The system falls back to default extensions from the repository's `config.json` (if any are defined). This is useful for development where default extensions may be enabled in the repository's config.
  - If `MDV_USER_CONFIG_PATH` is **not set**: Extensions are read from the repository's `config.json` file (if any are defined). This allows development environments to use default extensions without needing to set environment variables.
  
  This allows you to configure extensions at deployment time without modifying the repository configuration files or rebuilding the image. For development environments, you can add default extensions to the repository's `config.json` file, and they will be used automatically when `MDV_USER_CONFIG_PATH` is not set, or as a fallback when it is set but the specified file is missing or invalid. To use custom extensions in production, set `MDV_USER_CONFIG_PATH` and provide a valid user config file.

## Testing

### Running Tests

To run the standard test suite (excluding performance tests):

```bash
make test
```

To run only performance/stress tests:

```bash
make test-performance
```

To run all tests including performance tests:

```bash
poetry run pytest mdvtools/
```

### Performance Testing in CI

Performance tests are resource-intensive and are not run by default in CI to conserve resources. They can be triggered in several ways:

1. **Commit message trigger**: Include `[perf-test]` in your commit message
2. **Pull request title**: Include `[perf-test]` in the PR title
3. **Manual workflow dispatch**: Go to the GitHub Actions tab and manually trigger the "Python Performance Tests" workflow

Example commit message:
```
feat: add new conversion feature [perf-test]
```

The performance tests will run on large datasets (up to 100k cells × 8k genes) and verify:
- Memory efficiency during conversion
- Conversion speed benchmarks
- Memory leak detection
- Chunked normalization performance

### Test Markers

Tests are categorized using pytest markers:
- `@pytest.mark.performance`: Performance/stress tests (excluded from regular CI)
- `@pytest.mark.slow`: Slow-running tests
- `@pytest.mark.integration`: Integration tests
