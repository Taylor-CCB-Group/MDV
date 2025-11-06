from __future__ import annotations
from typing import TYPE_CHECKING, Optional
import pandas as pd
import regex as re
from mdvtools.llm.templates import packages_functions
import json
import ast
import subprocess
import tempfile
import os

if TYPE_CHECKING:
    from mdvtools.mdvproject import MDVProject


def extract_code_from_response(response: str):
    """Extracts Python code from a markdown string response."""
    # Use a regex pattern to match content between triple backticks
    code_pattern = r"```python(.*?)```"
    #matches = re.findall(code_pattern, response, re.DOTALL)
    match = re.search(code_pattern, response, re.DOTALL)

    if match:
        # Extract the matched code and strip any leading/trailing whitespaces
       return f"{match.group(1).strip()}"
    # we should at least issue a warning here, no good will come of this...
    return ""
    # if not matches:
    #     return ""
    # return matches[-1].strip()


def extract_explanation_from_response(response: str):
    """Extracts the explanation text by removing code blocks from the response."""
    # Remove all code blocks (```python...``` or ```...```)
    explanation = re.sub(r"```.*?```", "", response, flags=re.DOTALL)
    return explanation.strip()


def prepare_code(result: str, data: Optional[str | pd.DataFrame], project: MDVProject, log = print, 
                 modify_existing_project=False, view_name="default"):
    """Given a response from the LLM, extract the code and post-process it, 
    attempting to ensure that 
    - parameters are appropriately ordered.
    - it will not run a server for the project
    - etc
    In the longer term, the intention is to move away from this approach, which is liable to break as the template changes,
    or as soon as it gets into the hands of actual users...
    """
    original_script = extract_code_from_response(result)
    log("# Extracted code from the response, without any modifications:")
    log(original_script)

    # log("# Apply the reorder transformation")
    # modified_script = reorder_parameters(original_script, data)
    lines = original_script.splitlines()

    # Find the starting line index
    start_index = next(
        (i for i, line in enumerate(lines) if not line.strip().startswith(('import', 'from'))), None
    )

    if start_index is not None:
        # Capture all lines starting from the first 'def'
        captured_lines = "\n".join(lines[start_index:])
    else:
        log("Pattern not found")
        # this seems like it may actually be an error - and maybe not recoverable at this point
        # I am adding this line to placate pyright which is complaining about `captured_lines` being possibly unbound
        captured_lines = "# WARNING:::: No code captured from the response when calling prepare_code().\n"

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
        #final_code = final_code.replace("project.add_datasource(datasource_name, data_frame_obs)", "# project.add_datasource(datasource_name, data_frame_obs)")
        #final_code = final_code.replace("project.add_datasource(datasource_name_2, data_frame_var)", "# project.add_datasource(datasource_name_2, data_frame_var)")

        #final_code = final_code.replace("project.add_rows_as_columns_link", "# project.add_rows_as_columns_link")
        #final_code = final_code.replace("project.add_rows_as_columns_subgroup", "# project.add_rows_as_columns_subgroup") 

        # make sure the final code does not contain p.static
        final_code = final_code.replace("project.convert_to_static_page", "# project.convert_to_static_page")
        # all lines that include `data_frame` can be somewhat safely removed with the current template
        #final_code = re.sub(r".*data_frame.*", "", final_code) # commenting this line, it is causing problems when a histogram is being generated.
        final_code = final_code.replace("delete_existing=True", "delete_existing=False")
        # final_code = final_code.replace("\"default\"", f"\"{view_name}\"") # "default" was also used e.g. for `brush = "default"`
        ## assumption of `view_name = "default"` being present in the code no longer holds - often it will include a nice view name
        ## but sometimes there might be a problem with it clashing with existing views, or with quotes in the view name...
        # final_code = final_code.replace("view_name = \"default\"", f"view_name = \"{view_name}\"")
        # so we have a sticking plaster solution for this...
        final_code = patch_viewname(final_code, project)
        
        # Lint the code with ruff
        final_code = _lint_code_with_ruff(final_code, log)
        
    try:
        compile(final_code, "<string>", "exec") # will raise an exception if there is a syntax error
    except Exception as e:
        log(f"Error in the final code: {e}")
        log(final_code)
        # let it return anyway...
    return final_code

def _lint_code_with_ruff(code: str, log=print):
    """Formats the given code using ruff."""
    temp_file_path = None
    try:
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.py', encoding='utf-8') as temp_file:
            temp_file.write(code)
            temp_file_path = temp_file.name

        log("# Running ruff linting and fixing...")
        # Use --exit-zero to avoid raising an error if there are unfixable lint issues.
        ruff_result = subprocess.run(
            ['ruff', 'check', temp_file_path, '--fix', '--exit-zero'],
            capture_output=True, text=True
        )

        if ruff_result.stderr:
            log(f"# Ruff stderr:\n{ruff_result.stderr}")

        with open(temp_file_path, 'r', encoding='utf-8') as temp_file:
            linted_code = temp_file.read()
        log("# Ruff pass complete.")
        return linted_code

    except FileNotFoundError as e:
        log(f"# Ruff pass failed: {e}")
        # continue with the un-linted code
        return code
    finally:
        if temp_file_path and os.path.exists(temp_file_path):
            os.remove(temp_file_path)

def patch_viewname(code: str, project: MDVProject):
    """Given a code string, replace the view_name with a unique name, 
    attempting to escape any quotes that might have been in the original.
    """
    # Error: 'MDVProject' object is not callable... not sure where or why.
    view_name = parse_view_name(code)
    if view_name is None:
        raise ValueError("No view name found in the code.")
    print(f'original view_name: {view_name}')
    escaped_view_name = json.dumps(view_name) # this should escape any quotes in the view_name
    # but it also adds quotes around the view_name, so we need to remove them...
    escaped_view_name = escaped_view_name[1:-1]
    existing_views = set(project.views) # sometimes existing views do not capture all the existing view - trying a set instead of list to test if that works better [k for k in project.views]
    
    if view_name not in existing_views:
        # just in case the view_name isn't a duplicate, but might have had quotes in it
        print(f'patched view_name: {escaped_view_name}')
        complete_view_name = f"view_name = \"{view_name}\""
        escaped_complete_view_name = f"view_name = \"{escaped_view_name}\""
        return code.replace(complete_view_name, escaped_complete_view_name)
        #return code.replace(view_name, escaped_view_name)
    n = 1
    new_view_name = f"{escaped_view_name} ({n})"
    while new_view_name in existing_views:
        n += 1
        new_view_name = f"{escaped_view_name} ({n})"
    print(f'patched view_name: {new_view_name}')
    complete_view_name = f"view_name = \"{view_name}\""
    new_complete_view_name = f"view_name = \"{new_view_name}\""
    return code.replace(complete_view_name, new_complete_view_name)
    #return code.replace(view_name, new_view_name)

def parse_view_name(code: str):
    """Given a code string, extract the view_name from it using AST.
    This is more robust than using regex.
    """
    try:
        tree = ast.parse(code)
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                # We are looking for a simple assignment, e.g. view_name = "..."
                if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id == 'view_name':
                    # The value should be a string literal
                    if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                        view_name = node.value.value
                        return view_name
    except SyntaxError:
        # print("Failed to parse code with AST. It might contain a syntax error.")
        pass

    return None
