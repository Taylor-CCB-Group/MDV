from typing import TypedDict, Optional, Literal

"""
This file is a work in progress to define the structure of the chart configuration,
maybe TypedDicts would be a good way to define types of charts and their properties.

nb we should consider how to provide defaults for these values, 
where the responsibility lies for setting them, and how to validate them.

Also, I'm not sure the best way to deal with `"type": "wgl_scatter_plot"` etc.
I'd like it to be a generic or similar, and to be able to infer types from the type field...
I *don't* want to have to specify the "type" string of the chart when creating it.
"""


class BaseChart(TypedDict, total=False):
    """ """

    # type: str  # we could have a union of all possible types...
    title: str
    legend: str
    param: (
        str | list[str]
    )  # we could represent this with a more semantic type - also, may be different for different types of charts
    size: tuple[int, int]
    """
    When using `absolute` positioning, this will be the top-left corner of the chart in pixels.
    Ignored when `gridstack` positioning is used.
    """
    position: tuple[int, int]
    id: Optional[
        str
    ]  # we could generate these later if they are not provided (also check for duplicates)


class ColorLegend(TypedDict):
    display: bool
    pos: tuple[int, int]


class ColorProperties(TypedDict):
    color_by: str
    color_legend: ColorLegend


class AxisProperties(TypedDict, total=False):
    label: str
    size: int
    textSize: int
    tickfont: int


class Axes(TypedDict):
    x: AxisProperties
    y: AxisProperties


"""

from typing import TypedDict

class TypedDictWithDefaults(type(TypedDict)):  # type: ignore
    # A type dictionary with defaults.
    def __call__(self, *args, **kwargs):
        defaults = dict(
            a = 0, # set your defaults here
            b = 1
        )
        defaults.update(kwargs)
        return super().__call__(**defaults)


class A(TypedDictWithDefaults, metaclass=TypedDictWithDefaults):
    a: int
    b: int

print(A(), A(a=5, b=6))
# {'a': 0, 'b': 1} {'a': 5, 'b': 6}

"""

class ScatterChart(BaseChart, ColorProperties):
    type: Literal["wgl_scatter_plot"]
    axis: Axes

# class StackedRowChart(BaseChart, ColorProperties):
#     axis: Axes

def StackedRowChart(**kwargs: BaseChart | ColorProperties):
    return {"type": "stacked_row_chart", **kwargs}
