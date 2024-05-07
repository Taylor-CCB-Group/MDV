# Processing, Clustering and Visualising in MDV 3k PBMCs

To setup running MDV locally please follow the instructions below:

Clone the repository
`git clone https://github.com/Taylor-CCB-Group/MDV.git`\
`cd MDV`\
`git checkout pjt-dev`

From the MDV folder:

1. Install front-end dependencies:
    `npm i`
(If you don't have node installed so that you can run the npm command please download it from here: https://nodejs.org/en/download)

2. Setup Python virtual environment and build the front-end that it will use. On Unix-like systems, there is an npm script that will do this automatically, provided that you have Python 3.12 installed and Poetry is available in your PATH:
    `npm run python-setup`

This is equvalent to running the following commands, so in case it doesn't work run these:
    `python -m venv venv`\
    `source venv/bin/activate`\
    `cd python`\
    `poetry install --with dev`\
    `npm run build-flask-vite`

NOTE:
In case python does not run because the default is python3, try again the above commands replacing python with python3.
On Windows systems the `source venv/bin/activate` command will not work, instead you need to run `venv/Scripts/activate.bat`

In case poetry does not work, or you do not want to use it, install mdvtools using pip by:
`pip install -e python`

After the setup is done, please follow the instructions below:

1. In the open terminal run the command below. If you do not have jupyter notebook installed install it.\
    `jupyter notebook`
2. Open and run the notebook