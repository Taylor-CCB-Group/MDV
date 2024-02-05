import pandas as pd
from mdvtools.mdvproject import MDVProject
import os
import json

'''
Test what happens with NaN, Infinity...
Previously, although there was code for filtering `na = na[~numpy.isnan(na)]`,
this didn't help if there also happened to be Infinity.

The `json.dumps` default behaviour of `allow_nan=True` would then cheerfully output non-compliant JSON,
which contrary to what the Python documentation states is not, as of this writing (2024-02-05), 
compatible with *any* standard JavaScript based decoders I tested (let alone 'most').

We now use `allow_nan=False`, so the user would be alerted at project creation time rather than runtime,
but also correct the error earlier in the process so that the particular case of `isinf` is handled.

If there are no valid numbers at all in a numeric column, then an exception will be thrown,
in this example we catch that and as of now, end up with a datasource with 0 columns.
'''

project_folder = os.path.expanduser('~/mdv/bad_data')
p = MDVProject(project_folder, delete_existing=True)
df_good = pd.DataFrame({'a': [42]})
df_mixed = pd.DataFrame({'a': [42, float('nan'), float('inf')]})
df_bad = pd.DataFrame({'a': [float('nan')]})
try:
    p.add_datasource('good data', df_good)
    assert(len(p.datasources[0]['columns']) == 1)
    print('good data added ok')
    p.add_datasource('mixed data', df_mixed)
    assert(len(p.datasources[1]['columns']) == 1)
    print('mixed data added ok, column metadata should not contain any NaN/Infinity etc:')
    print(f"{json.dumps(p.datasources[1]['columns'][0])}")
    # expect `ValueError('zero-size array to reduction operation minimum which has no identity')`
    p.add_datasource('bad data', df_bad)
    assert(False) # df_bad is sufficiently degenerate that we don't expect to get this far.
except Exception as e:
    print(f"There was an exception - was it expected? (yes) Does it help us understand how to fix the problem? What state is the project datasource metadata in after this?")
    print(e)
    print(json.dumps(p.datasources, indent=2))


p.serve(port=5055)