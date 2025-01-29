from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.mdvproject import MDVProject
import pytest
import os
import tempfile
import scanpy as sc
from pathlib import Path


def get_h5ad_files():
    """Helper function to get all h5ad files from test_datasets folder."""
    test_folder = Path("test_datasets")
    if not test_folder.exists():
        raise FileNotFoundError(f"Test datasets folder not found at: {test_folder}")

    h5ad_files = list(test_folder.glob("*.h5ad"))
    if not h5ad_files:
        raise FileNotFoundError("No .h5ad files found in test_datasets folder")

    return h5ad_files


@pytest.fixture(params=get_h5ad_files())
def h5ad_file(request):
    """
    Parametrized fixture that yields each .h5ad file from the test_datasets folder.
    """
    file_path = request.param

    if not file_path.exists():
        raise FileNotFoundError(f"The .h5ad file was not found at: {file_path}")

    yield file_path


@pytest.fixture
def test_project_dir():
    """Fixture to create a temporary directory for the MDV project."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir


def test_convert_scanpy_to_mdv(h5ad_file, test_project_dir):
    """Test the conversion from Scanpy AnnData to MDV format."""
    # Read the h5ad file
    adata = sc.read(h5ad_file)

    # Convert to MDV format
    mdv = convert_scanpy_to_mdv(test_project_dir, adata)

    # Basic assertions
    assert isinstance(mdv, MDVProject)

    # Check if datasources were created using get_datasource_names
    datasource_names = mdv.get_datasource_names()
    assert "cells" in datasource_names, "Cells datasource not found"
    assert "genes" in datasource_names, "Genes datasource not found"

    # Get the datasources using proper method
    cell_datasource = mdv.get_datasource_metadata("cells")
    gene_datasource = mdv.get_datasource_metadata("genes")

    # Verify cell table
    assert cell_datasource is not None
    assert gene_datasource is not None

    # Test with label parameter
    label = "test_"
    mdv_labeled = convert_scanpy_to_mdv(
        test_project_dir, adata, label=label, delete_existing=True
    )
    labeled_datasource_names = mdv_labeled.get_datasource_names()

    # Check labeled datasources
    assert f"{label}cells" in labeled_datasource_names
    assert f"{label}genes" in labeled_datasource_names


def test_convert_scanpy_to_mdv_views(h5ad_file, test_project_dir):
    """Test the view preservation and creation functionality."""
    adata = sc.read(h5ad_file)

    # Test with delete_existing=True
    mdv_new = convert_scanpy_to_mdv(test_project_dir, adata, delete_existing=True)
    assert "default" in mdv_new.views
    assert "cells" in mdv_new.views["default"]["initialCharts"]
    assert mdv_new.views["default"]["initialCharts"]["cells"] == []

    # Test duplicate datasource creation (should raise FileExistsError)
    with pytest.raises(
        FileExistsError,
        match="Attempt to create datasource 'cells' failed because it already exists",
    ):
        convert_scanpy_to_mdv(test_project_dir, adata, delete_existing=False)

    # Test with delete_existing=False on a fresh directory
    new_project_dir = tempfile.mkdtemp()
    try:
        mdv_preserve = convert_scanpy_to_mdv(
            new_project_dir, adata, delete_existing=False
        )
        for view_name, view_data in mdv_preserve.views.items():
            assert "initialCharts" in view_data
            assert "cells" in view_data["initialCharts"]
            assert "genes" in view_data["initialCharts"]
    finally:
        # Clean up the temporary directory
        import shutil

        if os.path.exists(new_project_dir):
            shutil.rmtree(new_project_dir)
