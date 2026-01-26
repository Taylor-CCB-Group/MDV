from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

from spatialdata.models import get_table_keys
if TYPE_CHECKING:
    from spatialdata import SpatialData

ElementType = Literal["images", "labels", "points", "shapes", "tables"]
el_types = ["images", "labels", "points", "shapes", "tables"]

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

def sdata_to_mermaid(sdata: "SpatialData", orientation: Literal["LR", "TD"] = "LR") -> str:
    """
    Convert a SpatialData object to a Mermaid diagram.

    Shows the relationship between coordinate systems, elements, and annotations.
    """
    lines = [f'flowchart {orientation}']
    
    # Collect all coordinate systems and their elements
    all_cs_data = []
    all_annotations = {}  # (table_name, annotated_element) -> set of coordinate systems
    
    for cs_name in sdata.coordinate_systems:
        cs_data = sdata.filter_by_coordinate_system(cs_name)
        cs_elements = {}
        cs_annotations = {}
        for element_type in el_types:
            if hasattr(cs_data, element_type):
                collection = getattr(cs_data, element_type)
                for element_name, element in collection.items():
                    if element_type not in cs_elements:
                        cs_elements[element_type] = []
                    cs_elements[element_type].append(element_name)
                    if element_type == "tables":
                        annotated, _element_description, _instance_key = get_table_keys(element)
                        annotated_array = annotated if isinstance(annotated, list) else [annotated]
                        annotated_array = [a for a in annotated_array if a is not None]
                        cs_annotations[element_name] = annotated_array
                        # Track annotations for de-duplication
                        for annotated_element in annotated_array:
                            key = (element_name, annotated_element)
                            if key not in all_annotations:
                                all_annotations[key] = set()
                            all_annotations[key].add(cs_name)
        
        cs = CoordinateSystem(name=cs_name, elements=cs_elements, annotations=cs_annotations)
        all_cs_data.append((cs_name, cs))
    
    # Generate diagram lines for coordinate systems
    for cs_name, cs in all_cs_data:
        lines.append(_cs_to_mermaid(cs))
    
    # Add de-duplicated annotations (each table->element annotation appears only once)
    for (table_name, annotated_element), cs_set in all_annotations.items():
        lines.append(f'  {table_name} -->|annotates|{annotated_element}')
    
    result = "\n".join(lines)
    return f"```mermaid\n{result}\n```"