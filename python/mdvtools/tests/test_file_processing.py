from pathlib import Path

import pytest

from mdvtools.file_processing import (
    _read_rows_without_blanks,
    _validate_tabular_structure,
    ValidationError,
)


def write_tabular_file(tmp_path: Path, name: str, content: str) -> Path:
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")
    return path


def test_read_rows_without_blanks_respects_max_rows(tmp_path: Path):
    path = write_tabular_file(
        tmp_path,
        "many_rows.csv",
        "col1,col2\n\n1,2\n3,4\n5,6\n",
    )

    rows = _read_rows_without_blanks(str(path), ",", max_rows=2)

    assert rows == [["col1", "col2"], ["1", "2"]]


def test_validate_tabular_structure_allows_single_column_csv(tmp_path: Path):
    path = write_tabular_file(
        tmp_path,
        "single_column.csv",
        "sample_id\nA\nB\n",
    )

    _validate_tabular_structure(str(path), "single_column.csv", ",")


def test_validate_tabular_structure_rejects_duplicate_headers(tmp_path: Path):
    path = write_tabular_file(
        tmp_path,
        "duplicate_header.csv",
        "name,name\nA,B\n",
    )

    with pytest.raises(ValidationError, match="duplicate column names"):
        _validate_tabular_structure(str(path), "duplicate_header.csv", ",")


def test_validate_tabular_structure_rejects_inconsistent_rows(tmp_path: Path):
    path = write_tabular_file(
        tmp_path,
        "bad.tsv",
        "name\tage\nalice\t30\nbob\t40\tunexpected\n",
    )

    with pytest.raises(ValidationError, match="inconsistent columns"):
        _validate_tabular_structure(str(path), "bad.tsv", "\t")


def test_validate_tabular_structure_rejects_csv_content_named_txt(tmp_path: Path):
    path = write_tabular_file(
        tmp_path,
        "wrong_delimiter.txt",
        "id,name,age\n1,Alice,34\n2,Bob,29\n",
    )

    with pytest.raises(ValidationError, match="looks comma-separated, not tab-separated"):
        _validate_tabular_structure(str(path), "wrong_delimiter.txt", "\t")


def test_validate_tabular_structure_rejects_tsv_content_named_csv(tmp_path: Path):
    path = write_tabular_file(
        tmp_path,
        "wrong_delimiter.csv",
        "id\tname\tage\n1\tAlice\t34\n2\tBob\t29\n",
    )

    with pytest.raises(ValidationError, match="looks tab-separated, not comma-separated"):
        _validate_tabular_structure(str(path), "wrong_delimiter.csv", ",")
