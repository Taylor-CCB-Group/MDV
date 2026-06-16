from typing import Any

from mdvtools.project_protocols import CreateProjectPromptProject
from mdvtools.mdvproject import MDVProject
from mdvtools.markdown_utils import create_project_markdown, create_column_markdown
from mdvtools.llm.dataset_scale import ProjectScale, assess_project_scale
from mdvtools.llm.datasource_roles import (
    infer_datasource_roles,
    format_feature_table_field_policy,
    format_marker_gene_scanpy_fallback_policy,
    format_marker_ranking_viz_policy,
    format_metadata_column_schema_policy,
    format_mdv_first_data_access_policy,
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
from mdvtools.llm.datasource_roles import (
    infer_datasource_roles,
    categorical_field_ids_from_metadata,
)
from mdvtools.llm.column_field_resolve import build_expression_wrapper_token

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



def get_createproject_prompt_RAG(
    project: CreateProjectPromptProject,
    path_to_data: str,
    datasource_name: str,
    final_answer: str,
    question: str,
    compact: bool = False,
    scale: ProjectScale | None = None,
) -> str:
    """
    Constructs a RAG prompt to guide LLM code generation for creating MDV plots.
    Handles both standard and gene-related queries.
    When ``compact=True``, omits long policy blocks for local models (Ollama).
    """
    if scale is None:
        scale = assess_project_scale(project, path_to_data)
    mdv_first_policy = format_mdv_first_data_access_policy(
        scale, path_to_data, compact=compact
    )
    # Build markdown context for the selected datasource; fall back to whole project if needed
    try:
        ds_meta = project.get_datasource_metadata(datasource_name)
        context_md = f"## **{datasource_name}:** ({ds_meta['size']} rows)\n\n" + create_column_markdown(ds_meta["columns"])
    except Exception:
        context_md = create_project_markdown(project)
    roles = infer_datasource_roles(project)
    if compact:
        expr_lines = ""
        if roles.expressions:
            expr_lines = "\n".join(
                f"- `{e.datasource_name}` (name_column=`{e.name_column}`, default_subgroup=`{e.subgroup_key}`)"
                for e in roles.expressions
            )
        else:
            expr_lines = "- (none detected)"
        return (
            """
    Project Data Context:

    """
            + context_md
            + """

    Datasource roles:
    - Observation/metadata datasource: `"""
            + roles.obs_datasource
            + """`
    - Expression datasources (feature tables):
"""
            + expr_lines
            + """

"""
            + mdv_first_policy
            + """

    Context: {context}

    The provided scripts demonstrate MDV chart construction patterns. Follow this workflow:
    1. Use MDVProject(project_path, delete_existing=False) when editing this project — never delete or recreate it.
    2. Follow the MDV-first data access rules above — do not load `.h5ad` or full datasources for chart-only views.
    3. Load existing tables with `project.get_datasource_as_dataframe(<ds>, columns=[field_ids...])` when a pandas preview is needed.
    4. Build charts with mdvtools chart classes; set `params` from field ids (not display names).
    5. Print bounded markdown previews before `project.set_view(...)` when showing tabular results.
    6. Use injected CHATMDV_* constants when present; do not call `project.get_datasource_roles()`.

    Agent-suggested columns and chart types:
    """
            + final_answer
            + """

    User question: """
            + question
            + final_answer
            + """

    Update these variables in generated code:
    - project_path = '"""
            + project.dir
            + """'
    - data_path = '"""
            + path_to_data
            + """'
    - datasource_name = '"""
            + datasource_name
            + """'
    - view_name = a descriptive string for the visualization.

    Chart type requirements (use agent-suggested types when valid):
    - Abundance Box plot: three categorical columns (repeat if fewer available).
    - Box plot / Violin plot: one categorical + one numerical column.
    - Density Scatter plot: two numerical + one categorical column.
    - Dot plot / Heatmap: one categorical + numerical columns.
    - Histogram: one numerical column.
    - Multiline chart: one numerical + one categorical column.
    - Pie Chart / Row Chart / Wordcloud: one categorical column.
    - Row summary box / Table Plot: any column(s).
    - Sankey / Stacked row chart: two categorical columns (repeat if only one).
    - Scatter plot (2D): two numerical + one color column.
    - Scatter plot (3D): three numerical + one color column.
    - Selection dialog plot: any column.
    - Text box plot: text only.

    Output format:
    - Return only one fenced ```python code block with the complete runnable script.
    - Do not add markdown narrative before or after the code block.
    - For chart requests, create the view(s) described above; for text/table-first requests use bounded print(...) calls.


"""
        )
    feature_field_policy = format_feature_table_field_policy(roles)
    marker_gene_policy = format_marker_gene_scanpy_fallback_policy(path_to_data, scale)
    table_chart_param_policy = format_obs_table_chart_param_policy()
    marker_ranking_viz_policy = format_marker_ranking_viz_policy()
    hybrid_routing_policy = format_scanpy_hybrid_routing_policy()
    viz_consistency_policy = format_visualization_consistency_policy()
    no_hallucination_policy = format_no_hallucination_chart_policy()
    metadata_column_policy = format_metadata_column_schema_policy()
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

    2. Data Access (MDV-first; Scanpy last resort):
"""+mdv_first_policy+"""
        - Observation/metadata datasource: `"""+roles.obs_datasource+"""`
        - Wrapper-capable expression datasources (feature tables):
"""+expr_lines+"""
        - Use `project.get_datasource_as_dataframe(<datasource_name>, columns=[...])` with only the Field IDs you need.
          Each entry may be a metadata **field id** (e.g. cluster column on `cells`) or a **chart FieldName wrapper** string
          for expression columns on the **observation (row) datasource** (same format as chart `params`, e.g.
          `rna_expr|GENE(rna_expr)|12`). Do not pass wrappers to the feature-table datasource dataframe.

        - Marker genes and missing columns (Scanpy last resort only):
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
        - Column metadata schema and multi-gene heatmaps (ChatMDV):
"""+metadata_column_policy+"""
        - Use a chart class (e.g., DotPlot, BoxPlot, SelectionDialogPlot) and set `params = [...]` using selected **field ids** (see Parameter Handling).
            - The suggested columns and chart type are given by """+final_answer+"""
        - Convert the chart to JSON using `convert_plot_to_json(plot)`
        - **Before** `project.set_view(...)`, print a concise preview of any table, aggregate, or subset the charts depend on using **GitHub-flavored markdown tables** so the chat UI can render them (for example `print(df_result.head(40).to_markdown(index=False))` or `print(grouped.head(20).to_markdown())`). Do **not** use `DataFrame.to_string()` for chat previews. Keep rows/columns bounded so the output stays readable; this stdout is shown in the chat window.
        - **Before** `project.set_view(...)`, ensure saved charts use the **same data pipeline** as those printed previews (see "Visualization vs analysis consistency" below).
        - Set the view using `project.set_view(view_name, view_object)`
        - IMPORTANT: Chart objects (including `TablePlot`) do not support row-subsetting methods like `set_row_indices(...)`
          or invented filter setters such as `set_background_filter(...)`. Only call `set_*` (and other public) methods
          that exist on the chart class in `mdvtools.charts.*` (preflight validates this before execution).
        - **Axis styling:** On `BoxPlot`, `ViolinPlot`, `ScatterPlot`, and `DotPlot`, use only
          `plot.set_axis_properties("x", {{"label": "...", "textSize": 13, "tickfont": 10}})` and the same for `"y"`.
          Never call `set_x_axis` or `set_y_axis` on those classes (valid only on `HistogramPlot`, `HeatmapPlot`,
          `MultiLinePlot`, `AbundanceBoxPlot`).
          On `ScatterPlot` / `BoxPlot` / `DensityScatterPlot`, use `set_filter(...)`, not `set_on_filter` (`ScatterPlot3D` only).
        - **Chat preview tables:** For `groupby` summaries before `set_view`, use simple aggregates such as
          `.median()`, `.describe()`, or tuple form `agg(col=("col", "median"))`. Do not mix named tuple aggregations with
          bare `lambda` functions in the same `groupby(...).agg(...)` call.
          To show a subset, either:
            (a) create a filtered datasource when building a **new** project from raw data, or use the **narrow** ChatMDV
                exception in section 3 (`chat_rank_genes_result` via `add_datasource` for long-format Scanpy marker tables
                when editing a project), or
            (b) include a `SelectionDialogPlot` so the user can filter interactively in the UI.
        - Expression on an embedding **within a group** (general pattern):
            - Build a **full-dataset** `ScatterPlot` with embedding field ids in `params` and color via
              `set_color_by` / `set_default_color` using a **wrapper token** from `build_expression_wrapper_token`
              (see section 5).
            - When the user asks to focus on a categorical group, add a `SelectionDialogPlot` whose `params` list the
              relevant categorical **Field IDs** from Project Data Context (section 8); do **not** pass row indices into
              the scatter chart.
            - For numeric summaries of a subset (means, counts, distributions in chat), filter with
              `project.get_datasource_as_dataframe(...)` and bounded `print(...to_markdown())`; keep the saved scatter
              on the full dataset unless the user explicitly needs a persisted filtered datasource.
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
              “default_subgroup=...” lines, or injected `CHATMDV_EXPR_SUBGROUP_KEY`). Do **not** copy tutorial examples that use `rna_expr` unless that key is listed.
            - `<feature>` is a value from the feature table's label column (`name_column` from the link, e.g. `gene_ids` or `name` — match the project)
            - `<index>` is the row index of that feature in the feature table (0-based)
        - **ChatMDV injected constants (always in generated scripts):** Use `CHATMDV_OBS_DATASOURCE`, `CHATMDV_EXPR_DATASOURCE`, `CHATMDV_EXPR_NAME_COLUMN`, `CHATMDV_EXPR_SUBGROUP_KEY`, `CHATMDV_EXPRESSIONS`, and `CHATMDV_CATEGORICAL_FIELD_IDS` when present. Do **not** call `project.get_datasource_roles()` (does not exist). Prefer these constants over re-calling `infer_datasource_roles(project)` when editing an existing project. For categoricals, use `CHATMDV_CATEGORICAL_FIELD_IDS` or `categorical_field_ids_from_metadata(...)` — never `col['dtype']`.
        - **RowsAsColumnsExpression fields:** Use `expr.datasource_name` and `expr.name_column` / `expr.subgroup_key`. Never use `expr.feature_table`, `expr.feature_datasource`, or similar invented attributes.
        - Build wrappers with `build_expression_wrapper_token(subgroup_key, feature, index)`; do not hand-build f-strings unless necessary.
        - Example (wrapper expression):
            ```python
            if CHATMDV_EXPR_DATASOURCE is None:
                raise RuntimeError("No rows-as-columns expression link; cannot build gene wrappers.")
            feature = "GENE_OR_PROTEIN_NAME"
            df_var = project.get_datasource_as_dataframe(
                CHATMDV_EXPR_DATASOURCE, columns=[CHATMDV_EXPR_NAME_COLUMN]
            )
            feature_index = df_var[CHATMDV_EXPR_NAME_COLUMN].astype(str).tolist().index(feature)
            wrapper = build_expression_wrapper_token(
                CHATMDV_EXPR_SUBGROUP_KEY, feature, feature_index
            )
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
        - Add `SelectionDialogPlot` when interactive filtering materially helps answer the question—for example,
          expression on an embedding **within** a categorical group, or when the user asks to explore subsets
          interactively (pair with a full-dataset scatter colored by expression; see section 4).
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
            - Table Plot: Requires any column(s). Use only when the user explicitly asks for a **table in the view**, to **browse rows**, or to **list** values—not for questions about **relationships**, **correlations**, or **distributions** between numeric columns (use scatter/density/box/violin instead).
            - Text box plot: Requires no columns, just text.
        - **Chart type selection (when multiple types are suggested):**
            - Questions about **relationship**, **correlation**, **association**, or **comparison between two numeric columns** → prefer **Scatter plot (2D)** or **Density Scatter plot** (with a categorical color column when available). Do **not** save a `TablePlot` unless the user explicitly requests a table chart or row browsing in the view.
            - **Table Plot** is for interactive row-level browsing in MDV when explicitly requested; tabular answers in chat use markdown `print(...)` previews (section 4), not a saved `TablePlot`, unless the user wants the table persisted in the project view.  
            - Violin plot: Requires only one categorical column and one numerical column.  
            - Wordcloud: Requires one categorical column.
    Output format:
    - Return only one fenced ```python code block with the complete runnable script.
    - Do not add markdown narrative, bullet lists, or explanations before or after the code block.
    - For chart-oriented requests, the script must create the chart view(s) described above.
    - For textual/table-first requests, put answers in bounded print(...) calls inside the script.



"""
    )
    return prompt_RAG
