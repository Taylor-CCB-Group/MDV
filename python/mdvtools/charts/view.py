import sys

# You may also pick one without version check, of course
if sys.version_info < (3, 11):
    from typing_extensions import TypedDict, Required
else:
    from typing import TypedDict, Required

from typing import Literal, Any
# from .chart_dicts import BaseChart


class Panel(TypedDict, total=False):
    """This describes a panel in a view."""

    layout: Literal["gridstack", "absolute"]
    panelWidth: int | float


class View(TypedDict, total=False):
    """This metadata of a view to be added to a project views.json file."""

    initialCharts: Required[dict[str, list[Any]]]  # todo
    dataSources: dict[str, Panel]
    links: dict[str, list[Any]]  # todo
