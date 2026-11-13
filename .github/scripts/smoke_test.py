"""Check that an installed mdvtools works. Run it with the interpreter that has mdvtools installed.

    python .github/scripts/smoke_test.py 1.4.0          # after `pip install mdvtools`
    python .github/scripts/smoke_test.py 1.4.0 --app    # after `pip install "mdvtools[app]"`
"""

import argparse
import importlib
import subprocess
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

import pandas

import mdvtools
from mdvtools.mdvproject import MDVProject

# modules guarded by require_extra("app", ...) that do no work when imported
APP_MODULES = ["mdvtools.dbutils.dbmodels", "mdvtools.auth.auth0_provider", "mdvtools.llm.langchain_mdv"]

parser = argparse.ArgumentParser()
parser.add_argument("version")
parser.add_argument("--app", action="store_true", help="check the app extra is installed")
args = parser.parse_args()

installed = version("mdvtools")
if installed != args.version:
    sys.exit(f"installed mdvtools {installed}, expected {args.version}")

static = Path(mdvtools.__file__).parent / "static"
if not any(static.rglob("*.js")):
    sys.exit(f"no frontend build in {static}; build it before `uv build`")

subprocess.run([Path(sys.executable).parent / "mdvtools", "--help"], check=True, capture_output=True)

with tempfile.TemporaryDirectory() as folder:
    project = MDVProject(folder, delete_existing=True)
    project.add_datasource("cells", pandas.DataFrame({"x": [1.0, 2.0, 3.0], "type": ["a", "b", "a"]}))
    assert [ds["name"] for ds in project.datasources] == ["cells"]

if args.app:
    for module in APP_MODULES:
        importlib.import_module(module)
else:
    try:
        importlib.import_module(APP_MODULES[0])
        sys.exit(f"{APP_MODULES[0]} imported without the app extra")
    except ModuleNotFoundError as error:
        if 'pip install "mdvtools[app]"' not in str(error):
            raise

print(f"mdvtools {installed} {'[app] ' if args.app else ''}smoke test passed")
