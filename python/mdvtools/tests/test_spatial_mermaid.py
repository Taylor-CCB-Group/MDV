"""Shared-structure mermaid diagrams for the conversion report."""

import dask.array as da
import numpy as np
from spatialdata import SpatialData
from spatialdata.models import Image2DModel

from mdvtools.spatial.mermaid import format_shared_structures, group_by_structure, sdata_to_mermaid


def _images(*names: str) -> SpatialData:
    pixels = da.from_array(np.zeros((1, 4, 4), dtype=np.uint8))
    return SpatialData(
        images={
            name: Image2DModel.parse(pixels, dims=("c", "y", "x"))
            for name in names
        }
    )


def test_same_structure_is_one_diagram_regardless_of_element_order():
    first = _images("morphology", "extra")
    second = _images("extra", "morphology")

    assert sdata_to_mermaid(first) == sdata_to_mermaid(second)
    groups = group_by_structure([("a.zarr", first), ("b.zarr", second)])

    assert len(groups) == 1
    names, diagram = groups[0]
    assert names == ["a.zarr", "b.zarr"]
    assert diagram.count("```mermaid") == 1
    assert "morphology" in diagram
    assert "extra" in diagram


def test_unannotated_table_makes_its_own_structure():
    from anndata import AnnData

    plain = _images("morphology")
    with_table = _images("morphology")
    with_table.tables["qc"] = AnnData(np.ones((2, 1), dtype=np.float32))
    groups = group_by_structure([("plain.zarr", plain), ("qc.zarr", with_table)])

    assert [names for names, _diagram in groups] == [["plain.zarr"], ["qc.zarr"]]
    assert "qc" in groups[1][1]


def test_different_structures_stay_in_separate_groups():
    shared = _images("morphology")
    other = _images("morphology", "transcripts")
    report = format_shared_structures(
        [("one.zarr", shared), ("two.zarr", shared), ("other.zarr", other)]
    )

    assert report.index("#### 2 objects") < report.index("#### 1 object")
    assert "- `one.zarr`" in report
    assert "- `two.zarr`" in report
    assert "- `other.zarr`" in report
    assert report.count("```mermaid") == 2
    assert "transcripts" in report.split("#### 1 object", 1)[1]
