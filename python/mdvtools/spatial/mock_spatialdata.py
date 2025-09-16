from mdvtools.mdvproject import MDVProject
import spatialdata as sd


def create_mock_spatialdata():
    mock_sdata = sd.datasets.blobs()    
    return mock_sdata
