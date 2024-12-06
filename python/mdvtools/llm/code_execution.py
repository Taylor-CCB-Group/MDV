import subprocess
import warnings
import tempfile
import concurrent.futures


plt_code = """import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
fig, ax = plt.subplots()

fruits = ['apple', 'blueberry', 'cherry', 'orange']
counts = [40, 100, 30, 55]
bar_labels = ['red', 'blue', '_red', 'orange']
bar_colors = ['tab:red', 'tab:blue', 'tab:red', 'tab:orange']

ax.bar(fruits, counts, label=bar_labels, color=bar_colors)

ax.set_ylabel('fruit supply')
ax.set_title('Fruit supply by kind and color')
ax.legend(title='Fruit color')

plt.show()"""

def execute_plt_test():
    with concurrent.futures.ThreadPoolExecutor() as executor:
        try:
            future = executor.submit(execute_code, plt_code)
            return future.result()
        except Exception as e:
            print(f"execute_plt_test caught error: {str(e)}")

def execute_code(final_code: str, open_code=False, log: callable = print):
    """
    Executes the code provided in the `final_code` string, returning `True` if the execution was successful,
    of `False` if an unhandled exception was raised or the execution resulted in an error code being returned.

    Current implementation writes the code to a temporary file and runs it using `subprocess`.
    This detail could change; this function is a wrapper around the current implementation.
    We may for example sandbox the code in a separate container, or use a different method of execution.

    Args:        
        final_code: The code to be executed
        open_code: For local development, if true will open the code in Visual Studio Code
        log: A callable that will be used to log messages
    """
    # exec(final_code) - I don't understand why this doesn't work
    # "name 'MDVProject' is not defined"
    with tempfile.NamedTemporaryFile(suffix=".py") as temp:
        log(f"Writing to temp file: {temp.name}")
        temp.write(final_code.encode())
        temp.flush()
        if open_code:
            subprocess.run(["code", temp.name])
        # try:
        #     subprocess.run(["python", temp.name])
        # except subprocess.CalledProcessError as e:
        #     log(f"CalledProcessError: {str(e)}")
        # except Exception as e:
        #     log(f"Error: {str(e)}")    
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            stdout, stderr = run_subprocess(["python", temp.name])
            if stderr:
                log(f"Standard Error: {stderr}")
            for warning in w:
                if issubclass(warning.category, UserWarning):
                    log(f"UserWarning encountered: {warning.message}")
            if stdout:
                log(f"Standard Output: {stdout}")
                return True
            else:
                log("--- No standard output - overall execution failed ---")
                log(final_code)
                return False


def run_subprocess(command: list[str]) -> tuple[str, str]:
    try:
        # Capture standard output and standard error
        result = subprocess.run(
            command,
            # shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            # text=True
        )
        return result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        # Handle subprocess errors
        print(f"Subprocess error: {e}")
        print(f"Standard Output: {e.stdout}")
        print(f"Standard Error: {e.stderr}")
        return None, e.stderr
    except Exception as e:
        # Handle other exceptions
        print(f"General error: {e}")
        return None, str(e)
