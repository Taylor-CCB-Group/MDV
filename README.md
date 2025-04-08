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

Run these comands from Powershell (Windows as administrator) or a Terminal (MaxOS/Linux) to clone the repository:

```
git clone https://github.com/Taylor-CCB-Group/MDV.git
cd MDV
```

Now you are in the MDV folder install front-end dependencies:

```
npm i
npm run build-flask-vite
```

Now install Python libraries 

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

You can now add Charts (scatterplots, rowcharts etc) and Views (contains the charts) to allow visualisation and querying of the data. A further tutorial on creating simple charts and tables is shown in [this tutorial](docs/tutorials/scanpy_python_tutorial_1.md).

