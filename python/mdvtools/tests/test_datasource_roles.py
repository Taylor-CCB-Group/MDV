from mdvtools.llm.datasource_roles import infer_datasource_roles


class FakeProject:
    def __init__(self):
        self._names = ["cells", "rna", "protein"]

    def get_datasource_names(self):
        return self._names

    def get_links(self, datasource, filter=None):
        assert datasource == "cells"
        assert filter == "rows_as_columns"
        return [
            {
                "datasource": "rna",
                "link": {
                    "rows_as_columns": {
                        "name_column": "name",
                        "name": "rna",
                        "subgroups": {
                            "rna_expr": {"name": "rna_expr", "label": "rna_expr", "type": "sparse"},
                            "rna_logged_counts": {"name": "rna_logged_counts", "label": "rna_logged_counts", "type": "sparse"},
                        },
                    }
                },
            },
            {
                "datasource": "protein",
                "link": {
                    "rows_as_columns": {
                        "name_column": "name",
                        "name": "protein",
                        "subgroups": {
                            "protein_expr": {"name": "protein_expr", "label": "protein_expr", "type": "sparse"},
                            "protein_clr": {"name": "protein_clr", "label": "protein_clr", "type": "dense"},
                        },
                    }
                },
            },
        ]


def test_infer_datasource_roles_prefers_cells_and_expr_subgroups():
    roles = infer_datasource_roles(FakeProject())
    assert roles.obs_datasource == "cells"
    assert [e.datasource_name for e in roles.expressions] == ["rna", "protein"]
    assert roles.expressions[0].name_column == "name"
    assert roles.expressions[0].subgroup_key == "rna_expr"
    assert roles.expressions[1].subgroup_key == "protein_expr"

