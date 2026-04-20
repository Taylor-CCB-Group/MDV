from mdvtools.mdvproject import MDVProject
from typing import Any
from mdvtools.markdown_utils import create_project_markdown, create_column_markdown
from mdvtools.llm.datasource_roles import (
    infer_datasource_roles,
    format_feature_table_field_policy,
    format_marker_gene_scanpy_fallback_policy,
    format_marker_ranking_viz_policy,
    format_no_hallucination_chart_policy,
    format_obs_table_chart_param_policy,
    format_scanpy_hybrid_routing_policy,
    format_visualization_consistency_policy,
)
prompt_data = """
Your task is to:  
1. Identify the type of data the user needs (e.g., categorical, numerical, etc.) by inspecting the DataFrames provided.
2. Use only the two DataFrames provided:
   - df1: observation/metadata table (e.g. cells)
   - df2: linked feature table for wrapper-based expression (e.g. rna/genes/protein). The gene/feature label column may be `name`, `gene_ids`, or another id — inspect `df2.columns`; do not assume `name` unless present.
3. Column selection logic:
   - For non-expression queries: select columns from df1 only. Inspect df1, using df1.columns
   - For expression-related queries (e.g., expression of a gene/protein, comparison of features, highest expressing features):
       a. Resolve the feature-label column from `df2.columns` (see Datasource roles / name_column when available). Use ONLY values from that column for feature identifiers — do NOT use unrelated metadata columns unless explicitly instructed.
       b. If a specific feature is mentioned by the user, check if it exists in that label column.
           - If it does not exist, assume the user provided a synonym and the label column may use different ids (e.g. Ensembl). Attempt mapping or lookup within `df2` only.
           - If a match is found, return that label value.
           - If no match is found, ignore the requested gene and select from available labels in `df2`.
       c. If no feature is mentioned, select one or more feature labels from `df2`.
       d. Only use the resolved label column from df2 for feature strings — do not invent column names.
4. Always return the list of required columns as a quoted comma-separated string, like:
   - "col1", "col2"
   - Or for expression-related: "col", "feature_name"   (make sure "col" is from df1) 
5. For expression-related queries:
   - Return both df1 columns and the selected feature name (from the resolved df2 label column).
   - Only return the name as a string (e.g., "gene_name")—do not wrap it.
6. NEVER create new DataFrames or modify existing ones.
7. Ensure that the selected columns match the visualization requirements:  
    - Abundance Box plot: Requires three categorical columns.  
      - If only one categorical variable is available, return it three times.  
      - If two are available, return one of them twice.  
    - Box plot: Requires only one categorical column and one numerical column.  
    - Density Scatter plot: Requires two numerical columns and one categorical column.  
    - Dot plot: Requires only one categorical column and any number of numerical columns.  
    - Heatmap: Requires only one categorical column and any number of numerical columns.  
    - Histogram: Requires one numerical column.  
    - Multiline chart: Requires one numerical column and one categorical column.  
    - Pie Chart: Requires one categorical column.  
    - Row Chart: Requires one categorical column.  
    - Row summary box: Requires any column(s).  
    - Sankey plot: Requires two categorical columns.  
      - If only one categorical variable is available, return it twice.  
    - Scatter plot (2D): Requires two numerical columns and one any column for color.
    - Scatter plot (3D): Requires three numerical columns and one any column for color.  
    - Selection dialog plot: Requires any column.  
    - Stacked row chart: Requires two categorical columns.  
      - If only one categorical variable is available, return it twice.  
    - Table Plot: Requires any column(s).  
    - Text box: Requires no columns, just text.  
    - Violin plot: Requires only one categorical column and one numerical column.  
    - Wordcloud: Requires one categorical column.  
8. Important: Clearly separate the selected columns with quotes and commas.
9. The column names are case sensitive therefore return them as they are defined in the dataframe.
10. Output format:
   - First line: The word "fields" following by quoted, comma-separated list of column names.
   - Second line: The word "charts" following by quoted, comma-separated list of suitable chart types for the selected columns.
11. NEVER explain your reasoning.
"""

packages_functions = """import os
import pandas as pd
try:
    import scanpy as sc
    HAS_SCANPY = True
except Exception:
    sc = None
    HAS_SCANPY = False
from mdvtools.mdvproject import MDVProject
from mdvtools.conversions import convert_scanpy_to_mdv
from mdvtools.charts.density_scatter_plot import DensityScatterPlot
from mdvtools.charts.heatmap_plot import HeatmapPlot
from mdvtools.charts.histogram_plot import HistogramPlot
from mdvtools.charts.dot_plot import DotPlot
from mdvtools.charts.box_plot import BoxPlot
from mdvtools.charts.scatter_plot_3D import ScatterPlot3D
from mdvtools.charts.row_chart import RowChart
from mdvtools.charts.scatter_plot import ScatterPlot
from mdvtools.charts.abundance_box_plot import AbundanceBoxPlot
from mdvtools.charts.stacked_row_plot import StackedRowChart
from mdvtools.charts.ring_chart import RingChart
from mdvtools.charts.pie_chart import PieChart
from mdvtools.charts.violin_plot import ViolinPlot
from mdvtools.charts.multi_line_plot import MultiLinePlot
from mdvtools.charts.table_plot import TablePlot
from mdvtools.charts.wordcloud_plot import WordcloudPlot
from mdvtools.charts.text_box_plot import TextBox
from mdvtools.charts.row_summary_box_plot import RowSummaryBox
from mdvtools.charts.selection_dialog_plot import SelectionDialogPlot
from mdvtools.charts.sankey_plot import SankeyPlot

import json
import numpy as np
import sys
"""
# def load_data(path):
#     #Load data from the specified CSV file.
#     return pd.read_csv(path, low_memory=False)

# def convert_plot_to_json(plot):
#     #Convert plot data to JSON format.
#     return json.loads(json.dumps(plot.plot_data, indent=2).replace("\\\\", ""))



def get_createproject_prompt_RAG(project: MDVProject, path_to_data: str, datasource_name: str, final_answer: str, question: str) -> str:
    """
    Constructs a RAG prompt to guide LLM code generation for creating MDV plots.
    Handles both standard and gene-related queries.
    """
    # Build markdown context for the selected datasource; fall back to whole project if needed
    try:
        ds_meta = project.get_datasource_metadata(datasource_name)
        context_md = f"## **{datasource_name}:** ({ds_meta['size']} rows)\n\n" + create_column_markdown(ds_meta["columns"])
    except Exception:
        context_md = create_project_markdown(project)
    roles = infer_datasource_roles(project)
    feature_field_policy = format_feature_table_field_policy(roles)
    marker_gene_policy = format_marker_gene_scanpy_fallback_policy(path_to_data)
    table_chart_param_policy = format_obs_table_chart_param_policy()
    marker_ranking_viz_policy = format_marker_ranking_viz_policy()
    hybrid_routing_policy = format_scanpy_hybrid_routing_policy()
    viz_consistency_policy = format_visualization_consistency_policy()
    no_hallucination_policy = format_no_hallucination_chart_policy()
    expr_lines = ""
    if roles.expressions:
        expr_lines = "\n".join(
            f"- `{e.datasource_name}` (name_column=`{e.name_column}`, default_subgroup=`{e.subgroup_key}`)"
            for e in roles.expressions
        )
    else:
        expr_lines = "- (none detected)"
    prompt_RAG = (
    """
    Project Data Context:

    """ + context_md + """

    Context: {context}

    The provided scripts demonstrate how to generate various data visualizations using the `mdvtools` library in Python.

    Each script follows this workflow:

    1. Setup:
        - Initialize an MDVProject instance using the method: MDVProject(project_path, delete_existing=True) when creating a new project.
        - When modifying an existing project (the default in this chat context), do NOT recreate or delete the project.

    2. Data Access (scanpy optional):
        - If HAS_SCANPY is True and `data_path` points to a valid .h5ad file, you MAY load AnnData:
            adata = sc.read_h5ad(data_path)
            data_frame_obs = adata.obs
            data_frame_var = adata.var.assign(name=adata.var_names.to_list())
        - Exception: for **marker genes / top genes per cluster** when MDV tables lack needed columns, if `data_path`
          is a `.h5ad` file, you SHOULD use `sc.read_h5ad` and Scanpy as described under "Marker genes" below.
        - OTHERWISE (no scanpy or no .h5ad):
            - DO NOT use scanpy.
            - Use the existing MDV datasources already registered in the project:
                - Observation/metadata datasource (df1): `"""+roles.obs_datasource+"""`
                - Wrapper-capable expression datasources (feature tables): 
"""+expr_lines+"""
            - Use `project.get_datasource_as_dataframe(<datasource_name>)` to load these tables, or
              `project.get_datasource_as_dataframe(<datasource_name>, columns=[...])` to load a subset. Each entry may be a metadata **field id** (e.g. cluster column on `cells`) or a **chart FieldName wrapper** string for expression columns on the **observation (row) datasource** (same format as chart `params`, e.g. `rna_expr|GENE(rna_expr)|12`). Do not pass wrappers to the feature-table datasource dataframe.

        - Marker genes and missing columns (ChatMDV):
"""+marker_gene_policy+"""

    3. Datasource Registration:
        - When modifying an existing project, do **not** call `project.add_datasource(...)` for arbitrary new tables **except**
          the **ChatMDV marker-result** case: persist a long-format Scanpy marker `DataFrame` with
          `project.add_datasource('chat_rank_genes_result', marker_df, replace_data=True, add_to_view=view_name)` as
          described under "Marker ranking vs DotPlot" (Parameter Handling). No other ad-hoc datasources.
        - When creating a **new** project from raw data (not the default chat flow), normal `add_datasource` rules apply.
        - Hybrid Scanpy routing contract (ChatMDV):
"""+hybrid_routing_policy+"""

    4. Plot Construction:
        - No hallucinated datasources or columns (all chart types):
"""+no_hallucination_policy+"""
        - Use a chart class (e.g., DotPlot, BoxPlot, SelectionDialogPlot) and set `params = [...]` using selected **field ids** (see Parameter Handling).
            - The suggested columns and chart type are given by """+final_answer+"""
        - Convert the chart to JSON using `convert_plot_to_json(plot)`
        - **Before** `project.set_view(...)`, print a concise text preview of any table, aggregate, or subset the charts depend on (for example `print(df_result.head(40).to_string())` or grouped top-N). Keep rows/columns bounded so the output stays readable; this stdout is shown in the chat window.
        - **Before** `project.set_view(...)`, ensure saved charts use the **same data pipeline** as those printed previews (see "Visualization vs analysis consistency" below).
        - Set the view using `project.set_view(view_name, view_object)`
        - IMPORTANT: Chart objects (including `TablePlot`) do not support row-subsetting methods like `set_row_indices(...)`.
          To show a subset, either:
            (a) create a filtered datasource when building a **new** project from raw data, or use the **narrow** ChatMDV
                exception in section 3 (`chat_rank_genes_result` via `add_datasource` for long-format Scanpy marker tables
                when editing a project), or
            (b) include a `SelectionDialogPlot` so the user can filter interactively in the UI.
        - Visualization vs analysis consistency (ChatMDV):
"""+viz_consistency_policy+"""

    5. Parameter Handling:
        - The string """+final_answer+""" guides which columns and chart types to use.
        - **Critical:** For cell-level (obs) columns, each string in `params` must be the datasource **Field ID** exactly as shown in the Project Data Context tables (the **Field ID** column), NOT the display "Column Name" alone. MDV matches charts to data by internal `field` keys; display names may differ from field ids.
        - **Feature / genes datasources only:** For charts whose `initialCharts` key is the **`genes`** (or other feature) datasource, use **Field ID** strings from the Project Data Context **for that datasource**—for example `gene_ids` if listed there—not guessed names like `name` unless listed. Do **not** assume `gene_ids` is a valid Field ID on **`cells`**.
        - If you copy from `data_frame_obs.columns`, prefer the MDV field id for that column from the context table when they differ.
        - Table chart / obs vs Scanpy output (ChatMDV):
"""+table_chart_param_policy+"""
        - For marker table views on `chat_rank_genes_result`, build `params` from persisted datasource metadata after
          `add_datasource(..., replace_data=True, ...)` (or from `marker_df.columns` when equivalent), then verify each
          param token exists in that datasource field set before `set_view(...)`. If tokens do not resolve, keep bounded
          `print(...)` output and skip saving a broken chart.
        - Marker ranking vs DotPlot / Heatmap (ChatMDV):
"""+marker_ranking_viz_policy+"""
        - Feature table field compatibility (ChatMDV):
"""+feature_field_policy+"""
        - For wrapper-based expression (rows-as-columns), DO NOT treat features as plain obs columns.
          Use the FieldName wrapper format:
            `"<subgroup_key>|<feature>(<subgroup_key>)|<index>"`
          where:
            - `<subgroup_key>` **must** match a key from the rows-as-columns link for this project (see **section 2**
              “default_subgroup=...” lines). Do **not** copy tutorial examples that use `rna_expr` unless that key is listed.
            - `<feature>` is a value from the feature table's label column (`name_column` from the link, e.g. `gene_ids` or `name` — match the project)
            - `<index>` is the row index of that feature in the feature table (0-based)
        - Example (wrapper expression):
            ```python
            subgroup_key = "gs"  # use the default_subgroup from section 2 / Datasource roles for this project
            feature = "GENE_OR_PROTEIN_NAME"
            feature_index = data_frame_var[name_column].tolist().index(feature)
            wrapper = f"{{subgroup_key}}|{{feature}}({{subgroup_key}})|{{feature_index}}"
            ```

    6. Gene-Related Queries:
        - Only perform wrapper expression if you have a usable feature table for the modality (e.g. `rna`, `protein`) with its `name_column` and at least one subgroup key.
        - If none are available, proceed without wrapper expression and prefer non-expression charts.
        - For marker / top-N gene requests, follow **section 2 "Marker genes and missing columns"** (cells vs `genes`, Scanpy when `.h5ad` available).
        - When combining marker-gene results with charts, follow **section 4 "Visualization vs analysis consistency"** so the view matches Scanpy/MDV logic.

    7. Chat-first textual/table outputs:
        - For requests that primarily ask for textual summaries, mappings, rankings, annotations, or table-like listings,
          prioritize chat output over plot creation.
        - **Precedence vs marker persistence:** For questions focused on **interpretation**, **cell-type prediction**, or
          **biology** without an explicit request to **keep tables or charts in the saved project view**, use **bounded
          `print(...)` and markdown explanation only**—do **not** call `project.add_datasource('chat_rank_genes_result', ...)`
          or add saved-view charts solely to duplicate chat output. Use that `add_datasource` path only when the user
          clearly needs the marker table **in the MDV view** (see section 3 and "Marker ranking vs DotPlot").
        - In these cases, provide the answer in:
            (a) the markdown explanation text, and
            (b) bounded script stdout via `print(...)` before any `project.set_view(...)`.
        - Do NOT create `TextBox` or `TablePlot` by default when the user intent is primarily textual/table output.

    8. Selection dialog usage:
        - Add `SelectionDialogPlot` only when interactive filtering materially helps answer the question
          (for example, when the user asks to explore subsets interactively).
        - Do not add a selection dialog unconditionally.

    9. Follow-up phrasing and intent routing:
        - Resolve follow-up phrases like “those genes”, “same”, and “visualize that” to the most recent marker/DE result in context.
        - If no prior marker list exists, compute a bounded top-N marker table first, then plot from that explicit gene list.
        - Route intent consistently:
            - Marker ranking/stat requests -> text/table-first (optional datasource persistence when view storage requested).
            - Text box requests -> datasource table/text chart on the marker datasource.
            - Heatmap/Dot/Bubble/Violin of marker genes -> expression visualization from valid expression fields/wrappers.
            - Cluster count/distribution -> aggregate from `cells.<cluster_key>` with optional ring/bar chart.
            - Cell-type prediction -> text-first with uncertainty notes; only persist when one-value-per-cell mapping is explicit.

    10. Your Task:
        - Interpret the user question and decide **whether** a saved-view chart is needed, and if so **which** chart type
          fits the question: """+question+final_answer+"""
        - Use **field ids** from the Project Data Context (Field ID column) in `params` as appropriate:
            - Wrap only expression feature names using the `<subgroup>|<feature>(<subgroup>)|<index>` form as shown in section 5.
            - For all non-gene columns, use the exact **Field ID** string. Field ids are case sensitive.
        - Use formatted f-strings for all dynamic strings.
        - Generate a valid Python script that implements the analysis. **If** a saved-view chart or table is appropriate,
          create it with the MDVProject framework; **otherwise** the script may be print/markdown-only with an empty or
          minimal `set_view` (see sections 4–7 and "Marker ranking vs DotPlot").
        - Update these variables with these values:
            - project_path = '"""+project.dir+"""'
            - data_path = '"""+path_to_data+"""'  # may be empty; if empty or HAS_SCANPY is False, DO NOT use scanpy except for marker-gene workflows when this path is a `.h5ad` file (see section 2).
            - view_name = a string, in double quotes, describing what is being visualized.
            - datasource_name = '"""+datasource_name+"""'
        - The possible charts are given by """+final_answer+""" and should follow the following visualisation guidelines for each type of chart:
            - Abundance Box plot: Requires three categorical columns.  
                - If only one categorical variable is available, use it three times.  
                - If two are available, use one of them twice.  
            - Box plot: Requires only one categorical column and one numerical column.  
            - Density Scatter plot: Requires two numerical columns and one categorical column.  
            - Dot plot: Requires only one categorical column and any number of numerical columns.  
            - Heatmap: Requires only one categorical column and any number of numerical columns.  
            - Histogram: Requires one numerical column.  
            - Multiline chart: Requires one numerical column and one categorical column.  
            - Pie Chart: Requires one categorical column.  
            - Row Chart: Requires one categorical column.  
            - Row summary box: Requires any column(s).  
            - Sankey plot: Requires two categorical columns.  
                - If only one categorical variable is available, use it twice.  
            - Scatter plot (2D): Requires two numerical columns and one any column for color.
            - Scatter plot (3D): Requires three numerical columns and one any column for color.  
            - Selection dialog plot: Requires any column.  
            - Stacked row chart: Requires two categorical columns.  
                - If only one categorical variable is available, use it twice.  
            - Table Plot: Requires any column(s).  
            - Text box plot: Requires no columns, just text.  
            - Violin plot: Requires only one categorical column and one numerical column.  
            - Wordcloud: Requires one categorical column.
    Output format: Return python code for the requested result.
    - For chart-oriented requests, return code that creates the chart view(s).
    - For textual/table-first requests, return code that computes and prints concise, bounded tabular/text output in chat.

    After generating the code, include a detailed explanation in your response (renderable by markdown renderer) that covers:
    1. Why this chart is the best way to answer the question
    2. What biological insights can be gained from this visualization
    3. What subsequent analysis tasks could be performed based on these results

    The explanation will be displayed in the chat window automatically.



"""
    )
    return prompt_RAG
