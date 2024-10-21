# Multi Dimensional Viewer


![logo](images/mdv_logo.png)

The Multi-Dimensional Viewer (MDV) is a web-based application designed to help users analyse, annotate, and share multi-dimensional data from various sources (e.g., biological, statistical, etc.). MDV supports fast, interactive analysis even with large datasets (up to 10 million data items), thanks to its use of web workers, shared array buffers, and native arrays.
&nbsp;

![summary](images/summary.png)

## Key Features

* Large assortment of interactive charts/widgets may be used and embedded in an analysis. These include:
    * Scatter plots (2D and 3D)
    * Box plots
    * Heat maps
    * Histograms
    * Pie charts
    * Violin plots
    * Annotation tools
    * Interactive spatial biology charts 

* Multiple views or pages of data may be created to tell a story with the data

* Charts/Widgets can pop out into separate windows to take advantage of multiple screens

* Multiple data sources (tables) and can loaded and links defined between them

* Data can added and/or modified by the user

* Diverse range of data sources (API calls, static files) can be used by implementing custom data loaders 

* Runs in a web browser (installation not required for uploading and viewing data)

### System Requirements

- **Browser**: A modern browser (e.g., Chrome, Firefox, Safari, or Edge)
- **Memory**: At least 4 GB of RAM, even for handling large datasets (~10 million items), due to the lazy loading of data in raw bytes

You can browse existing projects at the [MDV Website](https://mdv.molbiol.ox.ac.uk/)
You can register at the [MDV Registration page](https://mdv.molbiol.ox.ac.uk/user/register?next=https://mdv.molbiol.ox.ac.uk/) for your own projects. 

We are working hard on a new version that will make uploading and installing data much easier. Please contact [Us](mailto:stephen.taylor@well.ox.ac.uk) if you want to upload your own data sets or wish to be added to the mailing list for news. 

### For Development or Running a Local Version from the Repository:

If you plan to contribute or run the latest development version, additional tools are required:
- **Browser**: A modern browser (e.g., Chrome, Firefox, Safari, or Edge)
- **Git**: For cloning the repository and version control
- **Node.js & npm**: For managing JavaScript dependencies ([Installation Guide](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm))
- **Python**: Version 3.10 or higher (3.12 is the most thoroughly tested and supported)
- **Poetry** (optional, but recommended): For managing Python dependencies, especially if you are contributing to the Python codebase

## Running On Local Machine

If you’re working with large datasets or want more control over your projects, you can install MDV locally. MDV is written in JavaScript designed to be embedded in a web page (https://mdv.molbiol.ox.ac.uk/). However in the python directory of this repository, there are some python scripts to format data to a specific file structure and compiled JavaScript that can display that format. There is also a lightweight server that runs locally to display projects

### Installation

#### From a GitHub release version

Run these comands from Powershell (Windows) or a Terminal (MaxOS/Linux) to clone the repository:

```
git clone https://github.com/Taylor-CCB-Group/MDV.git
cd MDV
```

Now you are in the MDV folder install front-end dependencies:

```
npm i
npm run build-flask-vite
```

Now install Python libaries 

## MacOS / Linux

```
python -m venv venv
source venv/bin/activate
pip install -e python
```

## Windows

```
python -m venv venv
venv/Scripts/activate
pip install -e python
```

### Running a test project

This example will build and run a project based on the `pbmc3k_processed` dataset from `scanpy`:

```
python python/mdvtools/test_projects/scanpy_pbmc3k.py
```

Note: This will open a internet browser on your machine at the URL "localhost:5052". Depending on how your machine's firewall is configured you may see an error saying the port is blocked. You can either allow this port to be used or edit the last line of the script to use a port that is unused e.g.
```
p.serve(port=8080)  
```

When the above script is run and the page is loaded you will see something like this:

![image](https://github.com/user-attachments/assets/db22272e-37f4-497b-a914-ee415508ca45)

You can now add Charts (scatterplots, rowcharts etc) and Views (contains the charts) to allow visualisation and querying of the data.

### Displaying a subset of the data from the original MDV publication (Single cell spatial analysis reveals inflammatory foci of immature neutrophil and CD8 T cells in COVID-19 lungs)

Download the data.

https://zenodo.org/record/6513508/files/hyp_example_data.zip?download=1

Then cd to the python directory in mdv

```
cd path/to/mdv/python
```

Open a python shell
```
python
```

Create an MDV project from the downloaded folder and display it in a browser:

```python
from mdvtools.mdvproject import MDVProject
p = MDVProject("/path/to/hyp_example_data")
p.serve()
```

This will open a browser window at http://localhost:5000/ but you will need to go to
[http://127.0.0.1:5000](http://127.0.0.1:5000) to avoid permissions errors. Note that this will fail if the front-end code has not been built after checking out the repository - `npm run build-flask-vite` to update it, or use a release version which should have the necessary build output already present in `python/templates/static`.

## Running on a server

The default data storage is an hdf5 file which is a compressed files that allows random read/write access to the data. However, it cannot be accessed directly but requires some kind of wrapper e.g. h5py. Hence access via http calls directly is not possible and backend code is required to display an MDV project in a web page. However, it can be converted to simple continuous compressed blocks of data for each column and a json index . This allows direct access via an http request with a range header:-

```python
from mdvtools.mdvproject import MDVProject
p = MDVProject("/path/to/hyp_example_data")
p.convert_to_static_page("/path/to/myapp/")
```

The function also creates a simple home page for the project (index.html), which can be customized. The folder can then be put on the server and accessed via:-
```
https://myserver.com/path/to/myapp
```

## Development

* clone the repository
* npm install

### Git blame without formatting commits

The `.git-blame-ignore-revs` file is used to ignore commits when running `git blame`. This is useful for ignoring formatting commits when trying to find the original author of a line of code. To use it, run:

```bash
git config -e
```

And add the following lines:

```toml
[blame]
	ignoreRevsFile = ".git-blame-ignore-revs"
	markIgnoredLines = true
```

Mark ignored lines will prepend a `?` to the blame commit hashes for indirectly blamed lines.

### note - this documentation is somewhat deprecated - [for dev-branch](#dev-branch)

You can run a project in development mode, which allows you do debug the JavaScript code and make changes, which will be reflected in the browser.

To use a project that has been converted into a static webpage (`convert_to_static_page()`) just specify the location of the folder when you start the dev server.

```
webpack serve --env dev=folder:/path/to/myproject
```
Or you can use data that is being served from a project. Run `serve()` on a project, which by default will start a server on local host at port 5050, and then start the development server

```
webpack serve --env dev=http://127.0.0.1:5050
```
In both cases the server will be running at localhost:8080

### Building the App
```
npm run build-flask-vite
```

#### note - this documentation is somewhat deprecated - [for dev-branch](#dev-branch)

This will build JavaScript that is is suitable for use with the 'static' project format and the lightweight inbuilt server in the python module. It puts the JavaScript files and assets in python/mdvtools/static/js and python/mdvtools/static/img respectively. When a static project is created, these files are copied over to the project's folder so that it can run independently.

```
webpack --env build=production
```
The above will build the JavaScript code and exposes ChartManager and a few utility methods for loading data, but it is up to the developer to write methods to load/save data to the app. The js files are put in dist/basic and any assets in dist/basic/images.


There are a few options to customize the production build:-
* mode - by default the mode is production which minifies and optimizes the bundled code. If development is specified, more verbose, easier to debug code is generated
* outputpath - the default is dist/basic , but any folder can be specified
* assetpath - this is the path where the app will look for assets, by default it is ./, so the generated images folder needs to be placed in the same directory as the html page. It can be changed to a relative or absolute path e.g. /static/assets/js, in which case the images folder will need to be places the relevant location on the server
* nofont - if true then fontawesome will not be included in the build, it will need to be imported in the html file

```
webpack --env build=production mode=development \
              outputpath=/path/to/myfolder \
              assetpath =/static/assets \
              nofont=true \
```

## Dev branch

The 'dev' branch is currently being used for development. It is automatically deployed to https://mdv-dev.netlify.app/ when a commit is made to the branch.

This documentation, and some aspects of how things are arranged, should be considered a work in progress.

The index page can load data either served by the `mdvtools` module, or with static data from another server (specified with a `static` flag in the URL). In either case, a `dir` parameter is used to specify the location of the data.

For example, the following URL will load data from a static server in the current dev deployment:

[`https://mdv-dev.netlify.app/?static&dir=https://mdvstatic.netlify.app/ytrap2`](https://mdv-dev.netlify.app/?static&dir=https://mdvstatic.netlify.app/ytrap2)


This will look for files `datasources.json` etc. as generated by python at the URL `https://mdvstatic.netlify.app/ytrap2`. A server running on localhost can also be used to serve this static data, as long as it sets appropriate CORS headers. At no point will any data loaded into the system in this way be uploaded to any server.

To work on developing the JS client code in this branch, we are using Vite rather than WebPack, and the dev server can be run with `npm run dev`. The `index.html` file is intended to be usable without requiring a custom JS entrypoint for a given project, although if desired, other configurations options are available.

This server runs on `localhost:5170` and behaves similarly to the Netlify deployment, but with Hot-Module-Reloading of compatible parts of the client code, and the ability to proxy to a local `mdvtools.mdv_desktop` server, which is intended as a simple way to have multiple local projects running in the same development environment.

The same project running in a local dev-server can be accessed at `http://localhost:5170/?static&dir=https://mdvstatic.netlify.app/ytrap2`.

### `mdvtools` servers

With the `mdvtools` module installed in a Python environment, there are two servers that can be run to serve MDV projects locally. They can be used either as main entrypoints, invoked from the command line, or from within other scripts.

There is a helper npm script `npm run python-setup` that will create a virtual environment and install the `mdvtools` module. This doesn't currently work on Windows, but the `mdvtools` module can be installed manually with `pip install -e .` in the `python` directory.

#### `mdvproject`

```bash
(venv) $ python -m mdvtools.mdvproject /path/to/project
```

This will start a server that serves the project at `http://localhost:5050/`, and automatically open it in the default browser.

#### `mdv_desktop`

```bash
(venv) $ python -m mdvtools.mdv_desktop
```

Starts a server that serves projects from a `~/mdv/` directory. This is intended to be a simple way to run multiple projects in the same development environment without needing to restart etc, and is not intended for production use. It will also read a `config.json` file in the `mdv` directory which currently only contains a list of locations for other projects outside of the `mdv` directory - for instance, on an external volume.

This is currently hard-coded to run on port `5051`, and can be used with or without the Vite dev-server (running without the dev-server requires an initial `npm run build-flask-vite` to build the JS used in Flask templates).

### Vite configuration

The Vite configuration is found in `vite.config.mts`, and will behave differently depending on environment variables.

#### `build` environment variable: `'production' | 'dev_pt' | 'desktop' | 'desktop_pt'`

As of this writing, the ones ending `_pt` generally relate to builds using newer dev-branch features, while `production` and `desktop` mirror legacy WebPack configurations.

