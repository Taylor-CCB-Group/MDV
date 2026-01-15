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

* Multiple data sources (tables) and can loaded and links defined between them, for instance refering to gene expression matrix data in cells

* Data can added and/or modified by the user (for example, to edit annotations)

* Diverse range of data sources (API calls, static files) can be used by implementing custom data loaders 

* Runs in a web browser (installation not required for uploading and viewing data)

### System Requirements

- **Browser**: A modern browser (e.g., Chrome, Firefox, Safari, or Edge)
- **Memory**: At least 4 GB of RAM, even for handling large datasets (~10 million items), due to the lazy loading of data in raw bytes

For information about using MDV and uploading your data please go to the [MDV Website](https://mdv.ndm.ox.ac.uk/).

To install MDV using Docker see [MDV Installation](https://mdv.ndm.ox.ac.uk/docs/installation/installation-manual) or the [local installation manual](docs/installation/installation-manual.md).

You can browse existing projects an old [MDV Projects Website](https://mdv.molbiol.ox.ac.uk/projects)

### For Development or Running a Local Version from the Repository:

The canonical way of developing the application as of this writing is to use a [devcontainer](https://containers.dev/), which standardises the development environment and dependencies.

In this case, the prerequisites are:
- **Git** for cloning the repository and version control
- **Docker** for containerising the application and related dependencies
- **VSCode or similar, with devcontainer extension** this will automatically detect the presence of devcontainer configuration and use it to configure the appropriate Python interpreter and other tooling, such as linting for front-end and back-end code.
- **WSL if running in Windows** there can be filesystem issues with Docker in Windows, and it is recommended to clone the repo into a WSL distro in this case.

```sh
git clone https://github.com/Taylor-CCB-Group/MDV.git
cd MDV
code .
```

This should display a prompt to re-open in a devcontainer, which will initially take a few minutes to build and then open a server on http://localhost:5055.

For front-end development, a vite dev-server can be started with `npm run dev` from a terminal in the container. The main front-page for that is currently accessed at http://localhost:5170/catalog_dev


### Alternative steps for running locally

If you prefer to manage tools differently - for a lighter-weight environment, non VSCode-based editors, or integration into other workflows etc, this is also possible. Further documentation can be found on the [MDV documentation website](https://mdv.ndm.ox.ac.uk/)
