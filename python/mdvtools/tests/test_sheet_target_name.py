"""Display-name choice when the cleaning sheet repeats a name."""

from mdvtools.examples.rename_softdelete_mdv_project_fields import (
    duplicated_name_keys,
    sheet_target_name,
)


def test_blank_names_are_not_duplicates():
    assert duplicated_name_keys(["", "Cluster", " cluster ", ""]) == {"cluster"}


def test_unique_final_name_wins_over_a_duplicated_new_name():
    target, rule = sheet_target_name(
        "UMAP 1",
        "UMAP Dimension 1",
        "rna:X_umap_mindist_0.25_1",
        duplicate_finals=set(),
        duplicate_news={"umap dimension 1"},
    )
    assert target == "UMAP 1"
    assert rule == "workbook"


def test_duplicated_final_uses_a_unique_new_name():
    target, rule = sheet_target_name(
        "Cluster",
        "Leiden resolution 0.6",
        "leiden_res0.6",
        duplicate_finals={"cluster"},
        duplicate_news=set(),
    )
    assert target == "Leiden resolution 0.6"
    assert rule == "workbook"


def test_duplicated_final_and_new_name_use_displayed_name():
    target, rule = sheet_target_name(
        "Cluster",
        "Cluster",
        "leiden_res0.6",
        duplicate_finals={"cluster"},
        duplicate_news={"cluster"},
    )
    assert target == "leiden_res0.6"
    assert rule == "displayed_name"


def test_blank_displayed_name_is_reported_when_both_sheet_names_are_duplicated():
    target, rule = sheet_target_name(
        "Cluster",
        "Cluster",
        "",
        duplicate_finals={"cluster"},
        duplicate_news={"cluster"},
    )
    assert target == ""
    assert rule == "displayed_name_blank"
