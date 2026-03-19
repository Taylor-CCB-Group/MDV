
This folder contains a .toml file to build a scaled down version of MDVTools which has the compiled version of the JavaScript frontend and a lightweight server for local viewing/editing. Only the necessary python modules and dependencies are included and it is pip installable

To build , first of all you need to build the javascript. In the main directory:-
```bash
npm install
npm run build-flask-vite
```
Make sure you have hatch installed:-
```bash
pip install hatch
```

Then cd to this directory (mdvtools_lite_build) and build:-
```bash
hatch build
```

To test, run the following script, which will create a virtual environment, activate it,  and install mdvtools in it.
```bash
python -m venv /path/to/venv
source /path/to/venv/bin/activate
pip install /path/to/mdv/mdvtools_lite_build
```

In the environment you can than run various tests e.g. create and view a spatial project
```bash
mdvtools convert-spatial --serve  path/to/outputproject  /path/to/spatialdata
```

To upload to pypi/testpypi make sure you have twine installed:-
```bash
pip install twine
```

and then upload (assuming you are in the mdvtools_lite_build directory)
```bash
twine upload --repository testpypi dist/* #upload to test pypi
twine upload  dist/* #upload to main repository
```

You will need a to be member of the mdvtools project and have the right permissions as you will need to authenticate

//TODO
include only necessary modules in the JavaScript build

