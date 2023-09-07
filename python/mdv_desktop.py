# from jupyterlab import labapp
from mdv.mdvproject import MDVProject

"""
for now, just a placeholder for something that will be an entry point for a 'desktop' version of MDV
i.e. something to allow users to easily explore their data via gui 
without having to install npm/python dependencies etc.
"""
# labapp.main()

MDVProject("/Users/petertodd/data/mdv_obj_test2").serve()