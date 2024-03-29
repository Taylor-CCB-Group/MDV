# Multi Dimensional Viewer


![logo](images/mdv_logo.png)

Multi Dimensional Viewer (MDV) is web based application for analyzing, annotating and sharing multi-dimensional data from different modalities.  It is inspired by [dc charts](https://dc-js.github.io/dc.js/) and [crossfilter](https://square.github.io/crossfilter/), but is performant with over 10 million data items due to the use of web workers, shared array buffers and native arrays.  
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

* Diverse range of data sources (API calls,static files) can be used by implementing custom data loaders 

* Runs in a web browser (installation not required for uploading and viewing data)


## Running On Local Machine

If you have large amounts of data or projects you may wish to install MDV locally. MDV is written JavaScript designed to be embedded in a web page (https://mdv.molbiol.ox.ac.uk/). However in the python directory of this repository, there are some python scripts to format data to a specific file structure and compiled JavaScript that can display that format. There is also a lightweight server that runs locally to display projects

### Installation

Download and unzip the repository

https://github.com/Taylor-CCB-Group/MDV/archive/refs/heads/main.zip

or clone it
```
git clone https://github.com/Taylor-CCB-Group/MDV.git
```


### System Requirements

* A modern browser
* python (3.6 or above)
* only 4GB of ram is required even for large datasets (~10 000 000 items) as data is lazily loaded as raw bytes

### Displaying example data
download the  data

https://zenodo.org/record/6513508/files/hyp_example_data.zip?download=1

Then cd to the python directory
```
cd path/to/mdv/python
```

Install `mdv` (using `editable` flag for development):

```
pip install -e .
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
[http://127.0.0.1:5000](http://127.0.0.1:5000) to avoid permissions errors.

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

You can run a project in development mode, which allows you do debug the JavaScript code and make changes, which will be reflected in the browser.

To use a project that has been converted into a static webpage (`convert_to_static_page()`) just specify the location of the folder when you start the dev server.

```
webpack serve --env dev=folder:/path/to/myproject
```
Or you can use data that is being served from a project. Run `serve()` on a project, which by default will start a server on local host at port 5000, and then start the development server

```
webpack serve --env dev=http://127.0.0.1:5000
```
In both cases the server will be running at localhost:8080

### Building the App
```
webpack --env build=desktop
```

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

The 'pjt-dev' branch is currently being used for development. It is automatically deployed to https://mdv-dev.netlify.app/ when a commit is made to the branch.

The index page is designed to load simple static data from another server, at a URL indicated by the `dir` parameter, for example:

`https://mdv-dev.netlify.app/?dir=https://mdvstatic.netlify.app/ytrap2`

This will look for files `datasources.json` etc. as generated by python at the URL `https://mdvstatic.netlify.app/ytrap2`. A server running on localhost can also be used to serve this static data, as long as it sets appropriate CORS headers. At no point will any data loaded into the system in this way be uploaded to any server.

There is an experimental Python script `csv_to_static.py` takes data in CSV format and generates the files needed for this static data mode. It currently assumes that there will also be a folder `images` - the operation of this script along with other associated ones may change in future and should not be considered stable in its current form.

To work on developing the JS client code in this branch, we are starting to use Vite rather than WebPack, and the dev server can be run with `npm run dev`. The `index.html` file is intended to be usable without requiring a custom JS entrypoint for a given project.
