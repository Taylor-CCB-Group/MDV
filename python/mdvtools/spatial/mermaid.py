from dataclasses import dataclass
from typing import TYPE_CHECKING

from spatialdata.models import get_table_keys
if TYPE_CHECKING:
    from typing import Literal
    from spatialdata import SpatialData

ElementType = Literal["images", "labels", "points", "shapes", "tables"]
el_types = ["images", "labels", "points", "shapes", "tables"]

@dataclass
class CoordinateSystem:
    name: str
    elements: dict[ElementType, list[str]]
    annotations: dict[str, list[str]] # table_name -> list of things it annotates

def cs_to_mermaid(cs: CoordinateSystem) -> str:
    lines = []
    for el_type, els in cs.elements.items():
        lines.append(f'  {cs.name} --> {cs.name}_{el_type}[{el_type}]')
        for el in els:
            lines.append(f'  {cs.name}_{el_type} -->{el}')
    for table_name, annotated_elements in cs.annotations.items():
        for annotated_element in annotated_elements:
            lines.append(f'  {table_name} -->|annotates|{annotated_element}')
    return "\n".join(lines)

def sdata_to_mermaid(sdata: SpatialData) -> str:
    """
    Convert a SpatialData object to a Mermaid diagram.

    Shows the relationship between coordinate systems, elements, and annotations.
    """
    lines = ['flowchart LR']
    for cs_name in sdata.coordinate_systems:
        cs_data  = sdata.filter_by_coordinate_system(cs_name)
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
                        annotated_arrary = annotated if isinstance(annotated, list) else [annotated]    
                        cs_annotations[element_name] = annotated_arrary
                    
     
        cs = CoordinateSystem(name=cs_name, elements=cs_elements, annotations=cs_annotations)
        lines.append(cs_to_mermaid(cs))
    return "\n".join(lines)