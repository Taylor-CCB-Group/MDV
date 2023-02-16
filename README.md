# Multi Dimensional Viewer


![logo](images/mdv_logo.png)

Multi Dimensional Viewer (MDV) is tool for analyzing, annotating  and sharing multi dimensional data.  It is inspired by [dc charts](https://dc-js.github.io/dc.js/) and [crossfilter](https://square.github.io/crossfilter/), but is performant with over 10 million data items due to the use of web workers, shared array buffers and native arrays.  
&nbsp;
![summary](images/summary.png)

## Key Features

* Large assortment of interactive charts/widgets such as spread sheets, genome browser,image viewer and 3d scatter plot

* Multiple views of the same data 

* Charts/Widgets can pop out into separate windows to take advantage of multiple screens

* Multiple data sources (tables) and can loaded and links defined between them

* Data can added and/or modified by the user

* Diverse range of data sources (API calls,static files) can be used by implementing custom data loaders 




## Running Locally

MDV is just JavaScript designed to be embedded in a web page (https://mdv.molbiol.ox.ac.uk/). However in the python directory of this repository, there are some python scripts to format data to a specific file structure and compiled JavaScript that can display that format. There is also a lightweight server that runs locally to display projects

### Installation

Download and unzip the repository

https://github.com/Taylor-CCB-Group/MDV/archive/refs/heads/main.zip

or clone it
```
git clone https://github.com/Taylor-CCB-Group/MDV.git
```


### System Requirements

* a modern browser
* python (3.6 or above)
* only 4GB of ram is required even for large datasets (~10 000 000 items) as data is lazily loaded as raw bytes

### Displaying example data
download the  data

https://zenodo.org/record/6513508/files/hyp_example_data.zip?download=1

Then cd to the python directory
```
cd path/to/mdv/python
```

Install the required python packages
```
pip install -r requirements.txt
```

Open a python shell
```
python
```

Create an MDV project from the downloaded folder and display it in a browser
```python
from mdv.mdvproject import MDVProject
p = MDVProject("/path/to/hyp_example_data")
p.serve()
```


## Development

* clone the repository
* npm install
* run the example - npm run run-ex --env example=basic_example

