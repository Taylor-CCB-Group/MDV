"""
Auto-generated Pydantic models from Zod schemas.
This file is generated automatically - do not edit manually.
Generated from TypeScript Zod schemas to ensure consistency.

To regenerate this file, run: npm run build-schemas
"""

from typing import List, Optional, Union, Dict, Any, Tuple
from pydantic import BaseModel, Field
from enum import Enum

# Data types that can be stored in columns
class DataType(str, Enum):
    INTEGER = "integer"
    DOUBLE = "double"
    TEXT = "text"
    TEXT16 = "text16"
    UNIQUE = "unique"
    MULTITEXT = "multitext"
    INT32 = "int32"

# Parameter types for chart configuration
class ParamType(str, Enum):
    TEXT = "text"
    NUMBER = "number"
    MULTITEXT = "multitext"
    TEXT16 = "text16"
    MULTI_COLUMN_NUMBER = "_multi_column:number"
    MULTI_COLUMN_ALL = "_multi_column:all"

# Serialized form of column queries
class RowsAsColsQuerySerialized(BaseModel):
    """Identifies this as a rows-as-columns query"""
    type: str = Field("RowsAsColsQuery", description="Identifies this as a rows-as-columns query")
    linkedDsName: str = Field(description="Name of the linked datasource that contains the row data")
    maxItems: int = Field(gt=0, description="Maximum number of items to retrieve from the linked datasource")

# Field specification can be a string (FieldName) or a serialized query
FieldSpec = Union[str, RowsAsColsQuerySerialized]

# GridStack layout configuration
class GridStackConfig(BaseModel):
    """GridStack layout configuration for dashboard positioning"""
    gssize: Optional[Tuple[int, int]] = Field(None, description="GridStack size as [width, height] in grid units")
    gsposition: Optional[Tuple[int, int]] = Field(None, description="GridStack position as [x, y] coordinates")

# Color configuration for columns
class ColorColumn(BaseModel):
    """Complex column object for color mapping"""
    field: str = Field(description="Column field identifier")
    name: str = Field(description="Human-readable column name")
    datatype: DataType = Field(description="Data type of the column")

# Color mapping configuration
class ColorBy(BaseModel):
    """Complex column object for color mapping"""
    column: ColorColumn = Field(description="Complex column object for color mapping")

# Base configuration that all charts extend
class BaseConfig(GridStackConfig):
    """Base configuration shared by all chart types, including layout, styling, and data mapping"""
    id: str = Field(description="Unique identifier for the chart instance")
    size: Tuple[float, float] = Field(description="Chart dimensions as [width, height] in pixels")
    title: str = Field(description="Display title for the chart")
    legend: str = Field(description="Legend text describing the chart's purpose")
    type: str = Field(description="Chart type identifier (e.g., 'scatter_plot', 'bar_chart')")
    param: List[FieldSpec] = Field(description="Array of field specifications defining the data columns used by this chart")
    title_color: Optional[str] = Field(None, description="CSS color value for the chart title")
    color_by: Optional[Union[str, ColorBy]] = Field(None, description="Field or column configuration used to determine color mapping")
    color_legend: Optional[Any] = Field(None, description="Custom color legend configuration")
    log_color_scale: Optional[bool] = Field(None, description="Whether to use logarithmic scaling for color values")
    trim_color_scale: Optional[Union[str, str]] = Field(None, description="Quantile trimming for color scale to handle outliers")
    color_overlay: Optional[float] = Field(None, description="Opacity value for color overlay effects")
    fallbackOnZero: Optional[bool] = Field(None, description="Whether to use fallback colors when values are zero")

# Chart-specific configuration schemas
class ScatterPlotConfig(BaseConfig):
    """Configuration for scatter plot charts showing relationships between two variables"""
    type: str = Field("scatter_plot", description="Scatter plot chart type")
    opacity: Optional[float] = Field(None, ge=0, le=1, description="Opacity of scatter plot points (0-1)")
    radius: Optional[float] = Field(None, gt=0, description="Radius of scatter plot points in pixels")

class BarChartConfig(BaseConfig):
    """Configuration for bar charts displaying categorical data"""
    type: str = Field("bar_chart", description="Bar chart type")

class HistogramConfig(BaseConfig):
    """Configuration for histogram charts showing data distribution"""
    type: str = Field("histogram", description="Histogram chart type")
    bins: Optional[int] = Field(None, gt=0, description="Number of bins for histogram data grouping")

class HeatmapConfig(BaseConfig):
    """Configuration for heatmap charts displaying matrix data with color intensity"""
    type: str = Field("heatmap", description="Heatmap chart type")

class DotPlotConfig(BaseConfig):
    """Configuration for dot plots showing individual data points"""
    type: str = Field("dot_plot", description="Dot plot chart type")

class BoxPlotConfig(BaseConfig):
    """Configuration for box plots showing data distribution statistics"""
    type: str = Field("box_plot", description="Box plot chart type")

class ViolinPlotConfig(BaseConfig):
    """Configuration for violin plots showing data density distribution"""
    type: str = Field("violin_plot", description="Violin plot chart type")

class PieChartConfig(BaseConfig):
    """Configuration for pie charts showing proportional data"""
    type: str = Field("pie_chart", description="Pie chart type")

class TableConfig(BaseConfig):
    """Configuration for data tables displaying tabular information"""
    type: str = Field("table", description="Table chart type")

class TextBoxConfig(BaseConfig):
    """Configuration for text boxes displaying static or dynamic text content"""
    type: str = Field("text_box", description="Text box chart type")
    text: str = Field(description="Text content to display in the text box")

class WordcloudConfig(BaseConfig):
    """Configuration for wordcloud charts showing text frequency"""
    type: str = Field("wordcloud", description="Wordcloud chart type")

class SankeyConfig(BaseConfig):
    """Configuration for Sankey diagrams showing flow between categories"""
    type: str = Field("sankey", description="Sankey diagram type")

class MultiLineChartConfig(BaseConfig):
    """Configuration for multi-line charts showing multiple time series or categories"""
    type: str = Field("multi_line_chart", description="Multi-line chart type")
    stacked: Optional[bool] = Field(None, description="Whether to stack multiple lines on top of each other")
    band_width: Optional[float] = Field(None, gt=0, description="Width of confidence bands around lines")
    intervals: Optional[int] = Field(None, gt=0, description="Number of intervals for data aggregation")
    scaletrim: Optional[bool] = Field(None, description="Whether to trim the scale to focus on data range")

class MultiBoxPlotConfig(BaseConfig):
    """Configuration for multi-box plots comparing distributions across categories"""
    type: str = Field("multi_box_plot", description="Multi-box plot chart type")

class AbundanceBoxPlotConfig(BaseConfig):
    """Configuration for abundance box plots showing distribution of abundance data"""
    type: str = Field("abundance_box_plot", description="Abundance box plot chart type")

class DensityScatterConfig(BaseConfig):
    """Configuration for density scatter plots showing point density with color intensity"""
    type: str = Field("density_scatter", description="Density scatter plot chart type")

class RowChartConfig(BaseConfig):
    """Configuration for row charts displaying data in horizontal bars"""
    type: str = Field("row_chart", description="Row chart type")

class StackedRowChartConfig(BaseConfig):
    """Configuration for stacked row charts with multiple data series per row"""
    type: str = Field("stacked_row_chart", description="Stacked row chart type")

class RowSummaryBoxConfig(BaseConfig):
    """Configuration for row summary boxes displaying aggregated statistics"""
    type: str = Field("row_summary_box", description="Row summary box chart type")

class SelectionDialogConfig(BaseConfig):
    """Configuration for selection dialogs allowing user interaction"""
    type: str = Field("selection_dialog", description="Selection dialog chart type")

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
    SelectionDialogConfig,
    BaseConfig
]

# View configuration
class ViewConfig(BaseModel):
    """Configuration for a complete view containing multiple datasources and their charts"""
    dataSources: Dict[str, Dict[str, Any]] = Field(description="Configuration for each datasource in the view")
    initialCharts: Dict[str, List[ChartConfig]] = Field(description="Initial charts to display for each datasource")

# Chart manager configuration
class ChartManagerConfig(BaseModel):
    """Top-level configuration for the chart manager controlling the entire dashboard"""
    initialCharts: Optional[List[ChartConfig]] = Field(None, description="Initial charts to create when the manager starts")
    all_views: Optional[List[str]] = Field(None, description="List of all available view names")
    current_view: Optional[str] = Field(None, description="Currently active view name")
    permission: Optional[str] = Field(None, description="Permission level for chart operations")
    gridstack: Optional[bool] = Field(None, description="Whether to use GridStack for chart layout")
    chat_enabled: Optional[bool] = Field(None, description="Whether chat functionality is enabled")
    mdv_api_root: Optional[str] = Field(None, description="Root URL for MDV API endpoints")
    onlyView: Optional[ViewConfig] = Field(None, description="Single view configuration when only one view is needed")

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
    return BaseConfig(**config_dict)

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
    
    chart_classes = {
        "scatter_plot": ScatterPlotConfig,
        "bar_chart": BarChartConfig,
        "histogram": HistogramConfig,
        "heatmap": HeatmapConfig,
        "dot_plot": DotPlotConfig,
        "box_plot": BoxPlotConfig,
        "violin_plot": ViolinPlotConfig,
        "pie_chart": PieChartConfig,
        "table": TableConfig,
        "text_box": TextBoxConfig,
        "wordcloud": WordcloudConfig,
        "sankey": SankeyConfig,
        "multi_line_chart": MultiLineChartConfig,
        "multi_box_plot": MultiBoxPlotConfig,
        "abundance_box_plot": AbundanceBoxPlotConfig,
        "density_scatter": DensityScatterConfig,
        "row_chart": RowChartConfig,
        "stacked_row_chart": StackedRowChartConfig,
        "row_summary_box": RowSummaryBoxConfig,
        "selection_dialog": SelectionDialogConfig,
    }
    
    chart_class = chart_classes.get(chart_type, BaseConfig)
    return chart_class(**base_config) 