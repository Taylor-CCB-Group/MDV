from mdvtools.llm.code_execution import execute_code

code_returns_error = """import sys
print("imported sys...")
sys.exit(1)
"""

print("testing execute_code with code that returns an error")
ok, stdout, stderr = execute_code(code_returns_error)
assert not ok, "expected error"

code_raises_exception = """raise Exception("test exception")"""
ok, stdout, stderr = execute_code(code_raises_exception)
assert not ok, "expected exception"

code_imports_mdv = """import mdvtools
print("imported mdvtools...")
"""
ok, stdout, stderr = execute_code(code_imports_mdv)
assert ok, "expect to be able to import mdvtools"

print('all `execute_code` tests passed') # pending more tests...