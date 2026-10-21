"""Check that an installed mdvtools works. Run it with the interpreter that has mdvtools installed.

    python .github/scripts/smoke_test.py 1.4.0
"""

import subprocess
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

import pandas

import mdvtools
from mdvtools.mdvproject import MDVProject

expected = sys.argv[1]
installed = version("mdvtools")
if installed != expected:
    sys.exit(f"installed mdvtools {installed}, expected {expected}")

static = Path(mdvtools.__file__).parent / "static"
if not any(static.rglob("*.js")):
    sys.exit(f"no frontend build in {static}; build it before `uv build`")

subprocess.run([Path(sys.executable).parent / "mdvtools", "--help"], check=True, capture_output=True)

with tempfile.TemporaryDirectory() as folder:
    project = MDVProject(folder, delete_existing=True)
    project.add_datasource("cells", pandas.DataFrame({"x": [1.0, 2.0, 3.0], "type": ["a", "b", "a"]}))
    assert [ds["name"] for ds in project.datasources] == ["cells"]

print(f"mdvtools {installed} smoke test passed")
