from mdvtools.llm.code_execution import execute_code

code_returns_error = """import sys
sys.exit(1)
"""

print("testing execute_code with code that returns an error")
ok = execute_code(code_returns_error)
assert not ok, "expected error"

code_raises_exception = """raise Exception("test exception")"""
ok = execute_code(code_raises_exception)
assert not ok, "expected exception"

