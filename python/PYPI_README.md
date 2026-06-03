# MDV Tools

**Multi-Dimensional Viewer (MDV)** is a web-based application for analysing and visualising
large, multi-dimensional data (single-cell, spatial, and general tabular data). `mdvtools`
is the Python package for **creating, manipulating, and viewing MDV projects** — including a
lightweight local server and the ability to export a project as a static website.

- Homepage / source: https://github.com/Taylor-CCB-Group/MDV
- Issues: https://github.com/Taylor-CCB-Group/MDV/issues

## Installation

```bash
pip install mdvtools
```

This installs the **slim core**: everything you need to build, manipulate, view, and serve
MDV projects, including spatial-data support.

For the **full application** — database-backed project management, the chat/LLM assistant,
and authentication — install the `app` extra:

```bash
pip install "mdvtools[app]"
```

If you call a full-app feature on a slim install, you'll get a clear message telling you to
`pip install "mdvtools[app]"` rather than a cryptic import error.

> **Testing a pre-release from TestPyPI:**
> ```bash
> pip install --index-url https://test.pypi.org/simple/ \
>             --extra-index-url https://pypi.org/simple/ mdvtools
> ```

## Quickstart

Create and serve an MDV project from the `pbmc3k` example dataset. `folder` is where the
project is stored — any location you can write to; it's created if it doesn't exist.

```python
import scanpy
from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.serverlite import serve_project

folder = "/path/to/mytestproject"
data = scanpy.datasets.pbmc3k_processed()

project = convert_scanpy_to_mdv(folder, data)
serve_project(project)            # serves at http://localhost:5050
```

Open http://localhost:5050 in your browser to explore the project. The project files live in
`folder` and can be re-served later or exported as a static site.

There's also a command-line interface:

```bash
mdvtools --help
```

## License

GPL-3.0-only. See the [repository](https://github.com/Taylor-CCB-Group/MDV) for full details.
