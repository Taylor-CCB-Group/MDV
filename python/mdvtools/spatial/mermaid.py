from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, Sequence

from spatialdata.models import get_table_keys
if TYPE_CHECKING:
    from anndata import AnnData
    from spatialdata import SpatialData

ElementType = Literal["images", "labels", "points", "shapes", "tables"]
el_types: list[ElementType] = ["images", "labels", "points", "shapes", "tables"]

@dataclass
class CoordinateSystem:
    name: str
    elements: dict[ElementType, list[str]]
    annotations: dict[str, list[str]] # table_name -> list of things it annotates

def _cs_to_mermaid(cs: CoordinateSystem) -> str:
    """
    Convert a coordinate system to mermaid diagram lines.
    
    Args:
        cs: internal CoordinateSystem object representation
    """
    cs_key = f"cs_{cs.name}"
    # add a declaration line which labels with the original name,
    # but with key being `cs_{name}` so that it isn't conflated with elements having same name
    lines = [f"  {cs_key}({cs.name})"]
    for el_type, els in cs.elements.items():
        lines.append(f'  {cs_key} --> {cs_key}_{el_type}[{el_type}]')
        for el in els:
            lines.append(f'  {cs_key}_{el_type} -->{el}')
    return "\n".join(lines)

def _coordinate_systems(sdata: "SpatialData") -> list[CoordinateSystem]:
    """Element names and table links in each coordinate system, in a stable order."""
    systems: list[CoordinateSystem] = []
    for cs_name in sorted(sdata.coordinate_systems):
        cs_data = sdata.filter_by_coordinate_system(cs_name)
        cs_elements: dict[ElementType, list[str]] = {}
        cs_annotations: dict[str, list[str]] = {}
        for element_type in el_types:
            collection = getattr(cs_data, element_type, None)
            if collection is None:
                continue
            for element_name in sorted(collection):
                cs_elements.setdefault(element_type, []).append(element_name)
                if element_type != "tables":
                    continue
                annotated, _element_description, _instance_key = get_table_keys(collection[element_name])
                annotated_array = annotated if isinstance(annotated, list) else [annotated]
                targets = sorted(str(target) for target in annotated_array if target is not None)
                cs_annotations[element_name] = targets
        systems.append(CoordinateSystem(name=cs_name, elements=cs_elements, annotations=cs_annotations))
    seen_tables = {name for cs in systems for name in cs.elements.get("tables", [])}
    missing_tables = sorted(name for name in sdata.tables if name not in seen_tables)
    if missing_tables:
        annotations: dict[str, list[str]] = {}
        for table_name in missing_tables:
            annotations[table_name] = _table_annotation_targets(sdata.tables[table_name])
        systems.append(
            CoordinateSystem(name="unassigned", elements={"tables": missing_tables}, annotations=annotations)
        )
    return systems


def _table_annotation_targets(table: "AnnData") -> list[str]:
    try:
        annotated, _region_key, _instance_key = get_table_keys(table)
    except ValueError:
        return []
    annotated_array = annotated if isinstance(annotated, list) else [annotated]
    return sorted(str(target) for target in annotated_array if target is not None)

def sdata_to_mermaid(sdata: "SpatialData", orientation: Literal["LR", "TD"] = "LR") -> str:
    """
    Convert a SpatialData object to a Mermaid diagram.

    Shows the relationship between coordinate systems, elements, and annotations.
    Element order is sorted, so two objects with the same structure produce the same diagram.
    """
    lines = [f'flowchart {orientation}']
    annotations: set[tuple[str, str]] = set()
    for cs in _coordinate_systems(sdata):
        lines.append(_cs_to_mermaid(cs))
        for table_name, targets in cs.annotations.items():
            annotations.update((table_name, target) for target in targets)

    for table_name, annotated_element in sorted(annotations):
        lines.append(f'  {table_name} -->|annotates|{annotated_element}')

    result = "\n".join(lines)
    return f"```mermaid\n{result}\n```"


def group_by_structure(
    named: Sequence[tuple[str, "SpatialData"]],
) -> list[tuple[list[str], str]]:
    """Group objects whose coordinate systems, elements, and table links match.

    Names stay in input order. A group is emitted when its first member appears.
    """
    groups: dict[str, list[str]] = {}
    order: list[str] = []
    for name, sdata in named:
        diagram = sdata_to_mermaid(sdata)
        if diagram not in groups:
            groups[diagram] = []
            order.append(diagram)
        groups[diagram].append(name)
    return [(groups[diagram], diagram) for diagram in order]


def format_shared_structures(named: Sequence[tuple[str, "SpatialData"]]) -> str:
    """One mermaid diagram per shared structure, listing the objects that use it."""
    blocks: list[str] = []
    for names, diagram in group_by_structure(named):
        noun = "object" if len(names) == 1 else "objects"
        listed = "\n".join(f"- `{name}`" for name in names)
        blocks.append(f"#### {len(names)} {noun}\n\n{listed}\n\n{diagram}")
    return "\n\n".join(blocks)
