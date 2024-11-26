import pandas as pd
import regex as re
from .templates import packages_functions

def extract_code_from_response(response: str):
    """Extracts Python code from a markdown string response."""
    # Use a regex pattern to match content between triple backticks
    code_pattern = r"```python(.*?)```"
    match = re.search(code_pattern, response, re.DOTALL)

    if match:
        # Extract the matched code and strip any leading/trailing whitespaces
        return match.group(1).strip()
    return None

def reorder_parameters(script: str, dataframe: str | pd.DataFrame):
    if isinstance(dataframe, str):
        df = pd.read_csv(dataframe)
    else:
        df = dataframe
    categorical_columns = df.select_dtypes(
        include=["object", "category"]
    ).columns.tolist()
    numerical_columns = df.select_dtypes(include=["number"]).columns.tolist()

    def is_categorical(column):
        return column in categorical_columns

    def is_numerical(column):
        return column in numerical_columns

    # Define a regex pattern to find function definitions that create BoxPlots
    # todo we could have an array of expected plot types and use it to build one big regex
    # that should make it easier to add new plot types.
    patterns = [
        re.compile(r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)BoxPlot\((.*?)\)", re.DOTALL),
        re.compile(r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)DotPlot\((.*?)\)", re.DOTALL),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)AbundanceBoxPlot\((.*?)\)",
            re.DOTALL,
        ),
        # re.compile(
        #     r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)HistogramPlot\((.*?)\)", re.DOTALL
        # ),
        re.compile(r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)RingChart\((.*?)\)", re.DOTALL),
        re.compile(r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)RowChart\((.*?)\)", re.DOTALL),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)StackedRowChart\((.*?)\)", re.DOTALL
        ),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)HeatmapPlot\((.*?)\)", re.DOTALL
        ),
    ]

    pattern_multiline = re.compile(r'def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)MultiLinePlot\((.*?)\)', re.DOTALL)

    for pattern in patterns:
        if pattern.search(script):
            # Define a regex pattern to find params and param patterns
            pattern_param = re.compile(r'params\s*=\s*\[.*?\]|param\s*=\s*".*?"')

            def reorder_params(match_param):
                matched_text = match_param.group(0)  # Get the entire matched text

                # Extract parameter names
                if "params" in matched_text:
                    param_list = re.findall(r"\'(.*?)\'", matched_text)
                    param_list = re.findall(r"\"(.*?)\"", matched_text)
                else:
                    param_list = [re.findall(r"\"(.*?)\"", matched_text)[0]]

                # Check for the presence of categorical and numerical variables
                has_categorical = any(is_categorical(param) for param in param_list)
                has_numerical = any(is_numerical(param) for param in param_list)

                # Add a categorical variable if none is present
                if not has_categorical and categorical_columns:
                    param_list.insert(0, categorical_columns[0])
                    has_categorical = True

                if len(param_list) < 2:
                    return matched_text  # No need to reorder if there are fewer than 2 parameters

                first_param = param_list[0]
                second_param = param_list[1]

                # Check the types of the parameters using the dataframe
                # if first_param in df.columns and second_param in df.columns:
                if has_categorical and has_numerical:
                    if not (is_categorical(first_param) and is_numerical(second_param)):
                        param_list[0], param_list[1] = param_list[1], param_list[0]

                # Reconstruct the parameters with reordered values
                if "params" in matched_text:
                    reordered_params = (
                        f"params = ['{param_list[0]}', '{param_list[1:]}']"
                    )
                else:
                    reordered_params = f'param = "{param_list[0]}"'

                return reordered_params.replace("'[", " ").replace("]'", "")

            # Substitute the matches with reordered parameters
            modified_script = re.sub(pattern_param, reorder_params, script)

            return modified_script

    if pattern_multiline.search(script):
        # Define a regex pattern to find params and param patterns
        pattern_param = re.compile(r'params\s*=\s*\[.*?\]|param\s*=\s*".*?"')
        
        def reorder_params_multiline(match_param):
            matched_text = match_param.group(0)  # Get the entire matched text

            # Extract parameter names
            if 'params' in matched_text:
                param_list = re.findall(r'\'(.*?)\'', matched_text)
                param_list = re.findall(r'\"(.*?)\"', matched_text)
            else:
                param_list = [re.findall(r'\"(.*?)\"', matched_text)[0]]
            
            # Check for the presence of categorical and numerical variables
            has_categorical = any(is_categorical(param) for param in param_list)
            has_numerical = any(is_numerical(param) for param in param_list)

            # Add a categorical variable if none is present
            if not has_categorical and categorical_columns:
                param_list.insert(0, categorical_columns[0])
                has_categorical = True

            if len(param_list) < 2:
                return matched_text  # No need to reorder if there are fewer than 2 parameters

            first_param = param_list[0]
            second_param = param_list[1]

            # Check the types of the parameters using the dataframe
            #if first_param in df.columns and second_param in df.columns:
            if has_categorical and has_numerical:
                if not (is_numerical(first_param) and is_categorical(second_param)):
                    param_list[0], param_list[1] = param_list[1], param_list[0]

            # Reconstruct the parameters with reordered values
            if 'params' in matched_text:
                reordered_params = f"params = ['{param_list[0]}', '{param_list[1:]}']"
            else:
                reordered_params = f'param = "{param_list[0]}"'

            return reordered_params.replace('\'[', ' ').replace(']\'','')

        # Substitute the matches with reordered parameters
        modified_script_multiline = re.sub(pattern_param, reorder_params_multiline, script)

        return modified_script_multiline

    return script

def prepare_code(result: str, data: str | pd.DataFrame, log: callable = print, modify_existing_project=False, view_name="default"):
    original_script = extract_code_from_response(result)

    log("# Apply the reorder transformation")
    modified_script = reorder_parameters(original_script, data)
    lines = modified_script.splitlines()

    # Find the starting line index
    start_index = next(
        (i for i, line in enumerate(lines) if not line.strip().startswith(('import', 'from'))), None
    )

    if start_index is not None:
        # Capture all lines starting from the first 'def'
        captured_lines = "\n".join(lines[start_index:])
    else:
        log("Pattern not found")

    # Log the prompt and the output of the LLM to the google sheets
    # log_to_google_sheet(sheet, str(context_information_metadata_name), output['query'], prompt_RAG, code)

    # logger('# Run the saved Python file. This will start a server on localhost:5050, open the browser and display the plot with the server continuing to run in the background.')
    # %run temp_code_3.pyc
    log("# Executing the code...")
    # - in order to get this to run in a chat context, we might want to get rid of the call to `p.add_datasource`
    final_code = f"""{packages_functions}\n{captured_lines}
else:
    main()"""
    final_code = final_code.replace("project.serve()", "# project.serve()")
    if modify_existing_project:
        # not at all robust... won't be needed in future
        final_code = final_code.replace("project.add_datasource", "# project.add_datasource")
        # make sure the final code does not contain p.static
        final_code = final_code.replace("project.convert_to_static_page", "# project.convert_to_static_page")
        # all lines that include `data_frame` can be somewhat safely removed with the current template
        final_code = re.sub(r".*data_frame.*", "", final_code)
        final_code = final_code.replace("delete_existing=True", "delete_existing=False")
        # final_code = final_code.replace("\"default\"", f"\"{view_name}\"") # "default" was also used e.g. for `brush = "default"`
        final_code = final_code.replace("view_name = \"default\"", f"view_name = \"{view_name}\"")
        
    return final_code
