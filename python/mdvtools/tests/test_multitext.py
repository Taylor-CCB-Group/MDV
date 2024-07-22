from mdvtools.mdvproject import MDVProject
from tempfile import TemporaryDirectory
import pandas as pd


def test_compliant_multitext():
    data = pd.DataFrame({'test_data': ['a', 'b', 'a; b', '']})

    # todo: TempProject more general pattern for tests like this / staging uploads?
    with TemporaryDirectory() as dir:
        p = MDVProject(dir)
        # as of this writing, an error is routinely raised internally when 't' is not already a datasource
        # p.add_datasource('t', data, columns=[{'name': 'test_data', 'datatype': 'text'}])
        p.add_datasource('t', data, columns=[{'name': 'test_data', 'datatype': 'text'}])
        different_cols = p.get_column('t', 'test_data') != data['test_data']
        assert(not isinstance(different_cols, bool))
        assert(not any(different_cols))
        # but the more serious issue is that it doesn't follow 'multitext' control flow properly,
        # and goes down a float32 path.
        p.add_datasource('mt', data, columns=[{'name': 'test_data', 'datatype': 'multitext', 'separator': ';'}])
        assert(p.datasources[1]['columns'][0]['datatype'] == 'multitext')
        different_cols = p.get_column("mt", "test_data") != data["test_data"]
        assert(not isinstance(different_cols, bool))
        assert(not any(different_cols))
        # this method is implicated in not gettting the right data type to send to the client
        # byte_data = p.get_byte_data([{ 'datasource': 'mt', 'column': 'test_data' }], None)  # <HDF5 dataset "__tags": shape (N,), type "<u2">
        # <u2 is the NumPy little-endian unsigned 2-byte integer datatype, seems right for text16
        # for 'text' <HDF5 dataset "__tags": shape (N,), type "|u1">
        # todo catch the error with __tags saved from frontend...


if __name__ == '__main__':
    test_compliant_multitext()