
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

To test, run the following script, which will create a virtual environment, activate it,  install mdvtools and serve an mdv project. Once you stop the server, the virtual environment will be deleted. You will need to change the paths to reflect your set-up.
```bash
python -m venv /path/to/venv
source /path/to/venv/bin/activate
pip install /path/to/mdv/mdvtools_lite_build
python -m mdvtools.serverlite /path/to/mdvproject
deactivate
rm -r /path/to/venv
```

//TODO
include only necessary modules in the JavaScript build

