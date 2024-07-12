import pandas as pd
import regex as re
import subprocess
import tempfile
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
    categorical_columns = df.select_dtypes(
        include=["object", "category"]
    ).columns.tolist()
    numerical_columns = df.select_dtypes(include=["number"]).columns.tolist()

    def is_categorical(column):
        return column in categorical_columns

    def is_numerical(column):
        return column in numerical_columns

    # Define a regex pattern to find function definitions that create BoxPlots
    patterns = [
        re.compile(r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)BoxPlot\((.*?)\)", re.DOTALL),
        re.compile(r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)DotPlot\((.*?)\)", re.DOTALL),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)AbundanceBoxPlot\((.*?)\)",
            re.DOTALL,
        ),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)HistogramPlot\((.*?)\)", re.DOTALL
        ),
        re.compile(r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)RingChart\((.*?)\)", re.DOTALL),
        re.compile(r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)RowChart\((.*?)\)", re.DOTALL),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)StackedRowChart\((.*?)\)", re.DOTALL
        ),
        re.compile(
            r"def\s+(\w*)\s*\((.*?)\):\s*\n(.*?)HeatmapPlot\((.*?)\)", re.DOTALL
        ),
    ]

    for pattern in patterns:
        if pattern.search(script):
            # Define a regex pattern to find params and param patterns
            pattern_param = re.compile(r'params\s*=\s*\[.*?\]|param\s*=\s*".*?"')

            def reorder_params(match_param):
                matched_text = match_param.group(0)  # Get the entire matched text

                # Extract parameter names
                if "params" in matched_text:
                    param_list = re.findall(r"\'(.*?)\'", matched_text)
                else:
                    param_list = [re.findall(r"\'(.*?)\'", matched_text)[0]]

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

    return script

def prepare_code(result: str, path_to_data: str, log: callable = print):
    original_script = extract_code_from_response(result)

    log("# Apply the reorder transformation")
    modified_script = reorder_parameters(original_script, path_to_data)
    lines = modified_script.splitlines()

    # Find the starting line index
    start_index = next(
        (i for i, line in enumerate(lines) if line.strip().startswith("def")), None
    )

    if start_index is not None:
        # Capture all lines starting from the first 'def'
        captured_lines = "\n".join(lines[start_index:])
    else:
        log("Pattern not found")

    # Log the prompt and the output of the LLM to the google sheets
    # log_to_google_sheet(sheet, str(context_information_metadata_name), output['query'], prompt_RAG, code)

    # logger('# Run the saved Python file. This will start a server on localhost:5050, open the browser and display the plot with the server continuing to run in the background.')
    # %run temp_code_3.py
    log("# Executing the code...")
    # - in order to get this to run in a chat context, we might want to get rid of the call to `p.add_datasource`
    final_code = f"""{packages_functions}\n{captured_lines}
else:
    main()"""
    final_code = final_code.replace("project.serve()", "# project.serve()")
    return final_code

def execute_code(final_code: str, open_code=False, log: callable = print):
    # exec(final_code) - I don't understand why this doesn't work
    # "name 'MDVProject' is not defined"
    with tempfile.NamedTemporaryFile(suffix=".py") as temp:
        log(f"Writing to temp file: {temp.name}")
        temp.write(final_code.encode())
        temp.flush()
        if open_code:
            subprocess.run(["code", temp.name])
        subprocess.run(["python", temp.name])
