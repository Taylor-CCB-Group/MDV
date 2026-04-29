"""
ChatMDV execution layer: runs generated Python in a subprocess.

This module is part of the ChatMDV boundary (see CHATMDV_BOUNDARY.md). It should surface
actionable diagnostics when subprocess execution fails, without changing core MDV data APIs.
"""
import subprocess
import warnings
import tempfile
import concurrent.futures
import queue
import threading
import time
from collections.abc import Callable


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


def _stderr_with_failed_source(stderr_diagnostic: str, final_code: str) -> str:
    """Append executed source so callers (e.g. chat error surfaces) can debug SyntaxError etc."""
    return (
        f"{stderr_diagnostic}\n\n"
        f"--- Failed Python source ({len(final_code)} chars) ---\n"
        f"{final_code}"
    )


def execute_plt_test():
    with concurrent.futures.ThreadPoolExecutor() as executor:
        try:
            future = executor.submit(execute_code, plt_code)
            return future.result()
        except Exception as e:
            print(f"execute_plt_test caught error: {str(e)}")

def execute_code(
    final_code: str,
    open_code=False,
    log=print,
    on_output_line: Callable[[str, str], None] | None = None,
    on_heartbeat: Callable[[float], None] | None = None,
):
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
            stdout, stderr = run_subprocess(
                ["python", temp.name],
                on_output_line=on_output_line,
                on_heartbeat=on_heartbeat,
            )
            if stderr:
                log(f"Standard Error: {stderr}")
            for warning in w:
                if issubclass(warning.category, UserWarning):
                    log(f"UserWarning encountered: {warning.message}")
            if stdout is not None:
                log(f"Standard Output: {stdout}")
                return True, stdout, stderr
            else:
                log("--- Execution failed ---")
                log(final_code)
                return False, stdout, _stderr_with_failed_source(stderr, final_code)


def run_subprocess(
    command: list[str],
    *,
    on_output_line: Callable[[str, str], None] | None = None,
    on_heartbeat: Callable[[float], None] | None = None,
    heartbeat_interval_s: float = 10.0,
    poll_interval_s: float = 0.25,
    time_fn: Callable[[], float] = time.monotonic,
) -> tuple[str | None, str]:
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        if process.stdout is None or process.stderr is None:
            return None, "Subprocess streams were not initialized."

        line_queue: queue.Queue[tuple[str, str | None]] = queue.Queue()
        stdout_chunks: list[str] = []
        stderr_chunks: list[str] = []

        def _pump_stream(source: str, stream) -> None:  # pragma: no cover - worker thread
            try:
                for line in iter(stream.readline, ""):
                    line_queue.put((source, line))
            finally:
                line_queue.put((source, None))

        stdout_thread = threading.Thread(
            target=_pump_stream, args=("stdout", process.stdout), daemon=True
        )
        stderr_thread = threading.Thread(
            target=_pump_stream, args=("stderr", process.stderr), daemon=True
        )
        stdout_thread.start()
        stderr_thread.start()

        closed_streams: set[str] = set()
        start = time_fn()
        last_activity = start
        last_heartbeat = start

        while len(closed_streams) < 2:
            try:
                source, line = line_queue.get(timeout=poll_interval_s)
                if line is None:
                    closed_streams.add(source)
                    continue
                last_activity = time_fn()
                if source == "stdout":
                    stdout_chunks.append(line)
                else:
                    stderr_chunks.append(line)
                if on_output_line:
                    on_output_line(source, line.rstrip("\n"))
            except queue.Empty:
                now = time_fn()
                if (
                    on_heartbeat is not None
                    and now - last_activity >= heartbeat_interval_s
                    and now - last_heartbeat >= heartbeat_interval_s
                ):
                    on_heartbeat(now - start)
                    last_heartbeat = now

        return_code = process.wait()
        stdout_text = "".join(stdout_chunks)
        stderr_text = "".join(stderr_chunks)
        if return_code == 0:
            return stdout_text, stderr_text

        diagnostic = (
            f"Subprocess failed with returncode={return_code}. "
            f"stdout_preview={stdout_text[:400]!r}. "
            f"stderr_preview={stderr_text[:400]!r}"
        )
        if stderr_text:
            diagnostic = f"{diagnostic}\n{stderr_text}"
        return None, diagnostic
    except subprocess.CalledProcessError as e:
        # Defensive fallback for callsites/tests that monkeypatch process runners with CPE.
        print(f"Subprocess error: {e}")
        print(f"Standard Output: {e.stdout}")
        print(f"Standard Error: {e.stderr}")
        decoded_out = e.stdout.decode(errors="replace") if isinstance(e.stdout, bytes) else (e.stdout or "")
        decoded_err = e.stderr.decode(errors="replace") if isinstance(e.stderr, bytes) else (e.stderr or "")
        diagnostic = (
            f"Subprocess failed with returncode={e.returncode}. "
            f"stdout_preview={decoded_out[:400]!r}. "
            f"stderr_preview={decoded_err[:400]!r}"
        )
        if decoded_err:
            diagnostic = f"{diagnostic}\n{decoded_err}"
        return None, diagnostic
    except Exception as e:
        # Handle other exceptions
        print(f"General error: {e}")
        return None, str(e)
