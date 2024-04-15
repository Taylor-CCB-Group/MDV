import pandas as pd
from mdvtools.mdvproject import MDVProject
import os
import json


def test_bad_data():
    """
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
    """

    path = os.path.join(os.path.dirname(__file__), "temp", "test_bad_data")
    if not os.path.exists(path):
        os.makedirs(path)
    p = MDVProject(path, delete_existing=True)
    df_good = pd.DataFrame({"a": [42]})
    df_mixed = pd.DataFrame({"a": [42, float("nan"), float("inf")]})
    # when explicitly setting a numeric column which we want to interpret as 'text' (so it can be a category, e.g. clusterID)
    df_bad_text = pd.DataFrame({"a": [0, 1, float("nan")]})
    df_bad = pd.DataFrame({"a": [float("nan")]})
    try:
        p.add_datasource("good data", df_good)
        assert len(p.datasources[0]["columns"]) == 1
        print("good data added ok")
        p.add_datasource("mixed data", df_mixed)
        assert len(p.datasources[1]["columns"]) == 1
        print(
            "mixed data added ok, column metadata should not contain any NaN/Infinity etc:"
        )
        print(f"{json.dumps(p.datasources[1]['columns'][0])}")

        # providing 'columns' without a 'datatype' key will cause a KeyError, and the column will be ignored
        # (although the datasource will still think it has a 'length' of 3 even though it has 0 columns)
        bad_columns = p.add_datasource(
            "bad column metadata", df_bad_text, columns=[{"name": "a", "type": "text"}]
        )
        assert len(bad_columns) == 1
        # this was co-ercing 'nan' to 0, which means that there is no difference between 'nan' and 0, which is bad.
        p.add_datasource(
            "bad text data", df_bad_text, columns=[{"name": "a", "datatype": "text"}]
        )
        assert (
            len(p.datasources[3]["columns"][0]["values"]) == 3
        )  # expect ['0', '1', 'nan']

        # expect `ValueError('zero-size array to reduction operation minimum which has no identity')`
        p.add_datasource("bad data", df_bad)
        assert (
            False
        )  # df_bad is sufficiently degenerate that we don't expect to get this far.
    except Exception as e:
        print(
            "There was an exception - was it expected? (yes) Does it help us understand how to fix the problem? What state is the project datasource metadata in after this?"
        )
        print(e)
        print(json.dumps(p.datasources, indent=2))

    # p.serve(port=5055)
