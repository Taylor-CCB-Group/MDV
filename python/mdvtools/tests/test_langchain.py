from mdvtools.llm.code_execution import execute_code


def test_success_with_output():
    """Tests that code that runs successfully and produces output is handled correctly."""
    print("testing execute_code with code that has output")
    code_with_output = """print("hello")"""
    ok, stdout, stderr = execute_code(code_with_output)
    assert ok, "expected success"
    assert stdout is not None
    assert stdout.strip() == "hello"
    assert stderr == "", f"stderr was not empty: {stderr}"


def test_success_no_output():
    """Tests that code that runs successfully and produces no output is handled correctly."""
    print("testing execute_code with code that has no output")
    code_no_output = """import os
x = 1"""
    ok, stdout, stderr = execute_code(code_no_output)
    assert ok, f"expected success even with no stdout, but got ok={ok}"
    assert stdout == ""
    assert stderr == "", f"stderr was not empty: {stderr}"


def test_returns_error():
    """Tests that code that exits with a non-zero status code is handled as a failure."""
    print("testing execute_code with code that returns an error")
    code_returns_error = """import sys
print("imported sys...")
sys.exit(1)"""
    ok, stdout, stderr = execute_code(code_returns_error)
    assert not ok, "expected error"
    # The current implementation of run_subprocess returns None for stdout on error.
    assert stdout is None
    # sys.exit(1) does not write to stderr.
    assert stderr == ""


def test_raises_exception():
    """Tests that code that raises an exception is handled as a failure."""
    print("testing execute_code with code that raises an exception")
    code_raises_exception = """raise Exception("test exception")"""
    ok, stdout, stderr = execute_code(code_raises_exception)
    assert not ok, "expected exception"
    assert stdout is None
    assert "Exception: test exception" in stderr
    assert "Traceback (most recent call last):" in stderr


def test_imports_mdv():
    """Tests that importing mdvtools works."""
    print("testing execute_code with mdvtools import")
    code_imports_mdv = """import mdvtools
print("imported mdvtools...")"""
    ok, stdout, stderr = execute_code(code_imports_mdv)
    assert ok, "expect to be able to import mdvtools"
    assert stdout is not None
    assert stdout.strip() == "imported mdvtools..."
    assert stderr == ""


def run_all_tests():
    """Runs all tests."""
    test_success_with_output()
    test_success_no_output()
    test_returns_error()
    test_raises_exception()
    test_imports_mdv()
    print("all `execute_code` tests passed")


if __name__ == "__main__":
    run_all_tests()