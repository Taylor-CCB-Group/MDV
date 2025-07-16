"""
Auto-generated Pydantic models from Zod schemas.
This file is generated automatically - do not edit manually.
Generated from TypeScript Zod schemas to ensure consistency.
"""

from typing import List, Optional, Union, Dict, Any, Tuple
from pydantic import BaseModel, Field
from enum import Enum

class datatype(str, Enum):
    INTEGER = "integer"
    DOUBLE = "double"
    TEXT = "text"
    TEXT16 = "text16"
    UNIQUE = "unique"
    MULTITEXT = "multitext"
    INT32 = "int32"


class ScatterPlotConfig(BaseModel):
    """Configuration for scatter plot charts showing relationships between two variables"""
    gssize: Optional[List[Any]] = Field(None, description="GridStack size as [width, height] in grid units")
    gsposition: Optional[List[Any]] = Field(None, description="GridStack position as [x, y] coordinates")
    id: str = Field(description="Unique identifier for the chart instance")
    size: List[Any] = Field(description="Chart dimensions as [width, height] in pixels")
    title: str = Field(description="Display title for the chart")
    legend: str = Field(description="Legend text describing the chart's purpose")
    type: str = Field(description="Scatter plot chart type")
    param: List[Union[str, Dict[str, Any]]] = Field(description="Array of field specifications defining the data columns used by this chart")
    title_color: Optional[str] = Field(None, description="CSS color value for the chart title")
    color_by: Optional[Any] = Field(None, description="Field or column configuration used to determine color mapping")
    color_legend: Optional[Dict[str, Any]] = Field(None, description="Custom color legend configuration")
    log_color_scale: Optional[bool] = Field(None, description="Whether to use logarithmic scaling for color values")
    trim_color_scale: Optional[str] = Field(None, description="Quantile trimming for color scale to handle outliers")
    color_overlay: Optional[float] = Field(None, description="Opacity value for color overlay effects")
    fallbackOnZero: Optional[bool] = Field(None, description="Whether to use fallback colors when values are zero")
    opacity: Optional[float] = Field(None, description="Opacity of scatter plot points (0-1)")
    radius: Optional[float] = Field(None, description="Radius of scatter plot points (units should be better defined)")


class BarChartConfig(BaseModel):
    """Configuration for bar charts displaying categorical data"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Bar chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class HistogramConfig(BaseModel):
    """Configuration for histogram charts showing data distribution"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Histogram chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")
    bins: Optional[float] = Field(None, description="Number of bins for histogram data grouping")


class HeatmapConfig(BaseModel):
    """Configuration for heatmap charts displaying matrix data with color intensity"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Heatmap chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class DotPlotConfig(BaseModel):
    """Configuration for dot plots showing individual data points"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Dot plot chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class BoxPlotConfig(BaseModel):
    """Configuration for box plots showing data distribution statistics"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Box plot chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")
    color_by: Optional[Any] = Field(None, description="")
    color_legend: Optional[Any] = Field(None, description="")
    log_color_scale: Optional[Any] = Field(None, description="")
    trim_color_scale: Optional[Any] = Field(None, description="")
    color_overlay: Optional[Any] = Field(None, description="")
    fallbackOnZero: Optional[Any] = Field(None, description="")


class ViolinPlotConfig(BaseModel):
    """Configuration for violin plots showing data density distribution"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Violin plot chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")
    color_by: Optional[Any] = Field(None, description="")
    color_legend: Optional[Any] = Field(None, description="")
    log_color_scale: Optional[Any] = Field(None, description="")
    trim_color_scale: Optional[Any] = Field(None, description="")
    color_overlay: Optional[Any] = Field(None, description="")
    fallbackOnZero: Optional[Any] = Field(None, description="")


class PieChartConfig(BaseModel):
    """Configuration for pie charts showing proportional data"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Pie chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class TableConfig(BaseModel):
    """Configuration for data tables displaying tabular information"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Table chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class TextBoxConfig(BaseModel):
    """Configuration for text boxes displaying static or dynamic text content"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Text box chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")
    text: str = Field(description="Text content to display in the text box")


class WordcloudConfig(BaseModel):
    """Configuration for wordcloud charts showing text frequency"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Wordcloud chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class SankeyConfig(BaseModel):
    """Configuration for Sankey diagrams showing flow between categories"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Sankey diagram type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class MultiLineChartConfig(BaseModel):
    """Configuration for multi-line charts showing multiple time series or categories"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Multi-line chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")
    stacked: Optional[bool] = Field(None, description="Whether to stack multiple lines on top of each other")
    band_width: Optional[float] = Field(None, description="Width of confidence bands around lines")
    intervals: Optional[float] = Field(None, description="Number of intervals for data aggregation")
    scaletrim: Optional[bool] = Field(None, description="Whether to trim the scale to focus on data range")


class MultiBoxPlotConfig(BaseModel):
    """Configuration for multi-box plots comparing distributions across categories"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Multi-box plot chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class AbundanceBoxPlotConfig(BaseModel):
    """Configuration for abundance box plots showing distribution of abundance data"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Abundance box plot chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class DensityScatterConfig(BaseModel):
    """Configuration for density scatter plots showing point density with color intensity"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Density scatter plot chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")
    color_by: Optional[Any] = Field(None, description="")
    color_legend: Optional[Any] = Field(None, description="")
    log_color_scale: Optional[Any] = Field(None, description="")
    trim_color_scale: Optional[Any] = Field(None, description="")
    color_overlay: Optional[Any] = Field(None, description="")
    fallbackOnZero: Optional[Any] = Field(None, description="")


class RowChartConfig(BaseModel):
    """Configuration for row charts displaying data in horizontal bars"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Row chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class StackedRowChartConfig(BaseModel):
    """Configuration for stacked row charts with multiple data series per row"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Stacked row chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class RowSummaryBoxConfig(BaseModel):
    """Configuration for row summary boxes displaying aggregated statistics"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Row summary box chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


class SelectionDialogConfig(BaseModel):
    """Configuration for selection dialogs allowing user interaction"""
    gssize: Optional[Any] = Field(None, description="")
    gsposition: Optional[Any] = Field(None, description="")
    id: Any = Field(description="")
    size: Any = Field(description="")
    title: Any = Field(description="")
    legend: Any = Field(description="")
    type: str = Field(description="Selection dialog chart type")
    param: Any = Field(description="")
    title_color: Optional[Any] = Field(None, description="")


# Union of all chart configuration types
ChartConfig = Union[
    ScatterPlotConfig,
    BarChartConfig,
    HistogramConfig,
    HeatmapConfig,
    DotPlotConfig,
    BoxPlotConfig,
    ViolinPlotConfig,
    PieChartConfig,
    TableConfig,
    TextBoxConfig,
    WordcloudConfig,
    SankeyConfig,
    MultiLineChartConfig,
    MultiBoxPlotConfig,
    AbundanceBoxPlotConfig,
    DensityScatterConfig,
    RowChartConfig,
    StackedRowChartConfig,
    RowSummaryBoxConfig,
    SelectionDialogConfig
]

class DataSource(BaseModel):
    """Complete configuration for a datasource including data structure, relationships, and visualization options"""
    name: str = Field(description="Human-readable name for the datasource")
    size: float = Field(description="Number of rows in the datasource")
    columns: List[Dict[str, Any]] = Field(description="Array of column configurations defining the data structure")
    columnGroups: Optional[List[Dict[str, Any]]] = Field(None, description="Logical groupings of related columns")
    links: Optional[Dict[str, Any]] = Field(None, description="Relationships with other datasources")
    images: Optional[Dict[str, Any]] = Field(None, description="Thumbnail images associated with datasource rows")
    large_images: Optional[Any] = Field(None, description="High-resolution images for detailed viewing")
    regions: Optional[Dict[str, Any]] = Field(None, description="Spatial regions configuration for spatial data")
    interactions: Optional[Any] = Field(None, description="Interaction data configuration")
    offsets: Optional[Dict[str, Any]] = Field(None, description="Spatial transformation configuration")
    genome_browser: Optional[Dict[str, Any]] = Field(None, description="Genome browser integration configuration")
    tree_diagram: Optional[Dict[str, Any]] = Field(None, description="Tree diagram visualization configuration")
    avivator: Optional[bool] = Field(None, description="Whether Avivator image viewer integration is enabled")
    row_data_loader: Optional[bool] = Field(None, description="Whether to use row-based data loading")
    binary_data_loader: Optional[bool] = Field(None, description="Whether to use binary data loading for performance")
    deeptools: Optional[Dict[str, Any]] = Field(None, description="DeepTools integration configuration")


class DataSourcesArray(BaseModel):
    """Array of multiple datasource configurations"""


# Structured output for data analysis
class DataAnalysisResult(BaseModel):
    """Result of structured data analysis"""
    selected_columns: List[str] = Field(description="List of selected column names from the dataframes")
    chart_types: List[str] = Field(description="List of suitable chart types for the selected columns")
    reasoning: str = Field(description="Brief explanation of the column selection")
    data_summary: Dict[str, Any] = Field(description="Summary statistics of the selected data")

# Structured output for chart generation
class ChartGenerationResult(BaseModel):
    """Result of structured chart generation"""
    chart_config: ChartConfig = Field(description="Generated chart configuration")
    python_code: str = Field(description="Generated Python code to create the chart")
    view_name: str = Field(description="Name for the generated view")
    dependencies: List[str] = Field(description="List of required Python packages")

# Complete structured output for a query
class StructuredQueryResult(BaseModel):
    """Complete structured result for a user query"""
    question: str = Field(description="Original user question")
    data_analysis: DataAnalysisResult = Field(description="Analysis of the data")
    chart_generation: ChartGenerationResult = Field(description="Generated chart")
    execution_result: Optional[Dict[str, Any]] = Field(None, description="Result of code execution")
    error: Optional[str] = Field(None, description="Error message if any")
    success: bool = Field(description="Whether the query was successful")

# Utility functions
def get_chart_type_from_config(config: ChartConfig) -> str:
    """Extract chart type from configuration"""
    if hasattr(config, 'type'):
        return config.type
    return "unknown"

def validate_chart_config(config_dict: Dict[str, Any]) -> ChartConfig:
    """Validate a chart configuration dictionary"""
    # This would need to be implemented based on the specific chart type
    # For now, we'll use the base config
    # ChartConfig is a Union, so we try each possible config class (there may be a better way)
    last_error = None
    for config_type in ChartConfig.__args__:
        try:
            return config_type(**config_dict)
        except Exception as e:
            last_error = e
    raise ValueError(f"Could not validate chart config: {last_error}")

def create_chart_config(
    chart_type: str,
    title: str,
    legend: str,
    params: List[str],
    size: Tuple[float, float] = (800, 600),
    **kwargs
) -> ChartConfig:
    """Create a chart configuration based on type and parameters"""
    base_config = {
        "id": f"chart_{hash(title)}",
        "title": title,
        "legend": legend,
        "param": params,
        "size": size,
        **kwargs
    }
    
    # ChartConfig is a Union, so we try each possible config class (there may be a better way)
    last_error = None
    for config_type in ChartConfig.__args__:
        try:
            return config_type(**base_config)
        except Exception as e:
            last_error = e
    raise ValueError(f"Could not create chart config: {last_error}")
