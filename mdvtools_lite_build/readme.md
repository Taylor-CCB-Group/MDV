
This folder contains a .toml file to build a scaled down version of MDVTools which has the compiled version of the JavaScript frontend and a lightweight server for local viewing/editing. Only the necessary python modules and dependencies are included and it is pip installable

To build , first of all you need to build the javascript. In the main directory:-
```
npm run vite-build
```
Make sure you have hatch installed:-
```
pip install hatch
```

Then cd to this directory (mdvtools_lite_build) and build:-
```
hatch build
```

To test run the following script, which will create a virtual environment, install in it the mdv version you have just built and create and serve a mdv project from the scanpy_pbm3k data. Once you stop the server, the virtual environment will be delete. You will need to change the paths to reflect your set up.
```bash
python -m venv /path/to/venv
source /path/to/venv/bin/activate
pip install /path/to/mdvtoolslite
python /path/to/mdv/python/test_projects/scanpy_pbm3k.py
deactivate
rm -r /path/to/venv
```

//TODO
include only necessary modules in the JavaScript build

