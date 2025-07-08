from mdvtools.llm.code_manipulation import parse_view_name


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