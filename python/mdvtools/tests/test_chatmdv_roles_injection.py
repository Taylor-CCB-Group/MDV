from mdvtools.llm.code_manipulation import prepare_code
from mdvtools.llm.code_preflight import validate_generated_code_preflight
from mdvtools.llm.datasource_roles import (
    build_chatmdv_roles_constants_block,
    collect_wrapper_subgroup_keys_for_project,
    infer_datasource_roles,
)
from mdvtools.mdvproject import MDVProject


class _FakeProject:
    def __init__(self):
        self.datasources = [{"name": "cells"}]
        self.dir = "/tmp/fake"
        self.views = []

    def get_datasource_names(self):
        return ["cells", "genes"]

    def get_links(self, datasource, filter=None):
        assert datasource == "cells"
        assert filter == "rows_as_columns"
        return [
            {
                "datasource": "genes",
                "link": {
                    "rows_as_columns": {
                        "name_column": "name",
                        "subgroups": {"gs": {"name": "gene_scores"}},
                    }
                },
            }
        ]

    def get_datasource_metadata(self, name):
        if name == "cells":
            return {
                "name": "cells",
                "columns": [
                    {"field": "leiden"},
                    {"field": "X_umap_1"},
                    {"field": "X_umap_2"},
                ],
                "links": {
                    "genes": {
                        "rows_as_columns": {
                            "subgroups": {"gs": {}},
                        }
                    }
                },
            }
        raise KeyError(name)


def test_build_chatmdv_roles_constants_block():
    block = build_chatmdv_roles_constants_block(_FakeProject())
    assert 'CHATMDV_OBS_DATASOURCE = "cells"' in block
    assert 'CHATMDV_EXPR_DATASOURCE = "genes"' in block
    assert 'CHATMDV_EXPR_NAME_COLUMN = "name"' in block
    assert 'CHATMDV_EXPR_SUBGROUP_KEY = "gs"' in block
    assert "CHATMDV_EXPRESSIONS" in block


def test_prepare_code_injects_roles_and_autofixes_bad_api():
    project = _FakeProject()
    llm_response = '''```python
def main():
    project = MDVProject("/tmp", delete_existing=False)
    roles = project.get_datasource_roles()
    expr = roles.preferred_expression()
    feature_table = expr.feature_table
```
'''
    code = prepare_code(
        llm_response,
        None,
        project,  # type: ignore[arg-type]
        log=lambda *_: None,
        modify_existing_project=True,
    )
    assert "CHATMDV_EXPR_DATASOURCE" in code
    assert "project.get_datasource_roles()" not in code
    assert "infer_datasource_roles(project)" in code
    assert ".feature_table" not in code
    assert ".datasource_name" in code


def test_preflight_rejects_hallucinated_roles_api():
    code = """
def main():
    project = MDVProject("/tmp", delete_existing=False)
    roles = project.get_datasource_roles()
"""
    res = validate_generated_code_preflight(code, allowed_wrapper_subgroup_keys={"gs"})
    assert res.ok is False
    assert any(i.code == "hallucinated_project_api" for i in res.issues)


def test_preflight_rejects_invalid_expression_attr():
    code = """
def main():
    expr = roles.preferred_expression()
    x = expr.feature_table
"""
    res = validate_generated_code_preflight(code)
    assert res.ok is False
    assert any(i.code == "invalid_expression_role_attr" for i in res.issues)


def test_collect_wrapper_subgroup_keys_pbmc3k():
    project = MDVProject("/app/mdv/pbmc3k_chat", delete_existing=False)
    keys = collect_wrapper_subgroup_keys_for_project(project)
    assert "gs" in keys
    roles = infer_datasource_roles(project)
    block = build_chatmdv_roles_constants_block(project)
    assert roles.preferred_expression() is not None
    assert "RALY" in block or "genes" in block
