from mdvtools.llm.code_manipulation import parse_view_name
from mdvtools.llm.code_manipulation import prepare_code, _defines_function_named_main


class TestParseViewName:
    def test_simple_assignment(self):
        code = 'view_name = "My cool view"'
        assert parse_view_name(code) == "My cool view"

    def test_assignment_with_comment(self):
        code = 'view_name = "Another view"  # this is a comment'
        assert parse_view_name(code) == "Another view"

    def test_assignment_with_extra_whitespace(self):
        code = 'view_name   =   "Spaced out view"'
        assert parse_view_name(code) == "Spaced out view"

    def test_no_assignment(self):
        code = 'project.add_view("a view")'
        assert parse_view_name(code) is None

    def test_syntax_error(self):
        code = 'view_name = "unterminated string'
        assert parse_view_name(code) is None

    def test_assignment_in_function(self):
        code = """
def create_view():
    view_name = "A view in a function"
    return view_name
"""
        assert parse_view_name(code) == "A view in a function"

    def test_different_quotes(self):
        code = "view_name = 'single quotes'"
        # The current implementation only supports double quotes because of how it reconstructs the string
        # for replacement, and how the regex was structured. The AST parser can handle it,
        # but the replacement logic in `patch_viewname` might be brittle.
        # For now, let's assert that it finds it, but a developer should be aware of this.
        # The prompt only specified finding the assignment, not fixing the whole chain.
        # Let's see if the current ast parser can handle single quotes.
        # The ast parser *should* handle it. Let's test it.
        # The value of an ast.Str node is just the string content.
        assert parse_view_name(code) == "single quotes"

    def test_multiline_code(self):
        code = """
import pandas as pd

def main():
    view_name = "a view name in multiline code"
    print(view_name)
"""
        assert parse_view_name(code) == "a view name in multiline code" 


def test_prepare_code_does_not_append_else_when_llm_includes_else_main(tmp_path):
    class FakeProject:
        def __init__(self):
            self._views = {}

        @property
        def views(self):
            return self._views

    llm = """```python
import os

def main():
    view_name = "v"
    return

if __name__ == "__main__":
    main()
else:
    main()
```"""
    out = prepare_code(llm, data=None, project=FakeProject(), modify_existing_project=False)
    # Must compile; regression for stray `else:` being appended.
    compile(out, "<string>", "exec")


def test_defines_function_named_main_helper():
    assert _defines_function_named_main("def main():\n    pass\n")
    assert _defines_function_named_main("async def main():\n    pass\n")
    assert _defines_function_named_main("def foo():\n    pass\ndef main():\n    pass\n")
    assert not _defines_function_named_main("x = 1\nprint(x)\n")
    assert not _defines_function_named_main("def run():\n    pass\n")


def test_prepare_code_appends_main_only_when_def_main_exists():
    class FakeProject:
        def __init__(self):
            self._views = {}

        @property
        def views(self):
            return self._views

    llm = """```python
def main():
    x = 1
```"""
    out = prepare_code(llm, data=None, project=FakeProject(), modify_existing_project=False)
    assert 'if __name__ == "__main__":\n    main()' in out


def test_prepare_code_multiline_add_datasource_compiles_with_modify_existing():
    """Regression: naive commenting of project.add_datasource broke multiline calls (unmatched ')')."""

    class FakeProject:
        def __init__(self):
            self._views = {}

        @property
        def views(self):
            return self._views

    llm = """```python
view_name = "v"
project.add_datasource(
    'chat_rank_genes_result',
    df,
    replace_data=True,
    add_to_view=view_name,
)
```"""
    out = prepare_code(llm, data=None, project=FakeProject(), modify_existing_project=True)
    compile(out, "<string>", "exec")
    assert "project.add_datasource(" in out
    assert "# project.add_datasource" not in out


def test_prepare_code_no_append_main_for_top_level_only():
    class FakeProject:
        def __init__(self):
            self._views = {}

        @property
        def views(self):
            return self._views

    llm = """```python
x = 1
print(x)
```"""
    out = prepare_code(llm, data=None, project=FakeProject(), modify_existing_project=False)
    assert 'if __name__ == "__main__":\n    main()' not in out
    compile(out, "<string>", "exec")