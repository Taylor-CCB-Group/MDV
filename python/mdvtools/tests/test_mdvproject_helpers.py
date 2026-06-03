import logging
from unittest.mock import MagicMock, patch

import pytest

from mdvtools import mdvproject as mp


def test_agent_debug_log_uses_logger_not_hardcoded_file():
    mock_logger = MagicMock(spec=logging.Logger)
    mp._agent_debug_log(
        "run-1",
        "H1",
        "test:loc",
        "hello",
        {"k": "v"},
        log=mock_logger,
    )
    mock_logger.debug.assert_called_once()
    mock_logger.warning.assert_not_called()


def test_agent_debug_log_logs_warning_on_logger_failure():
    mock_logger = MagicMock(spec=logging.Logger)
    mock_logger.debug.side_effect = RuntimeError("boom")
    mp._agent_debug_log("run-1", "H1", "test:loc", "hello", {}, log=mock_logger)
    mock_logger.warning.assert_called_once()


def test_resolve_rows_as_columns_prefers_candidate_over_other_link_single_subgroup():
    class Stub:
        def get_links(self, datasource, filter=None):
            del datasource
            return [
                {
                    "datasource": "genes",
                    "link": {
                        "rows_as_columns": {
                            "subgroups": {
                                "only": {"name": "wrong", "type": "dense"},
                            }
                        }
                    },
                },
                {
                    "datasource": "protein",
                    "link": {
                        "rows_as_columns": {
                            "subgroups": {
                                "rna_expr": {"name": "expr_mat", "type": "sparse"},
                            }
                        }
                    },
                },
            ]

    name, sparse = mp.MDVProject._resolve_rows_as_columns_subgroup(
        Stub(), "cells", "rna_expr"
    )
    assert name == "expr_mat"
    assert sparse is True


def test_resolve_rows_as_columns_single_subgroup_fallback_only_when_one_link():
    class OneLink:
        def get_links(self, datasource, filter=None):
            del datasource
            return [
                {
                    "datasource": "genes",
                    "link": {
                        "rows_as_columns": {
                            "subgroups": {
                                "gs": {"name": "scores", "type": "dense"},
                            }
                        }
                    },
                }
            ]

    name, sparse = mp.MDVProject._resolve_rows_as_columns_subgroup(
        OneLink(), "cells", "missing_key"
    )
    assert name == "scores"
    assert sparse is False


def test_resolve_rows_as_columns_no_fallback_when_one_link_has_many_subgroups():
    class OneLinkManySubgroups:
        def get_links(self, datasource, filter=None):
            del datasource
            return [
                {
                    "datasource": "genes",
                    "link": {
                        "rows_as_columns": {
                            "subgroups": {
                                "a": {"name": "mat_a", "type": "dense"},
                                "b": {"name": "mat_b", "type": "dense"},
                            }
                        }
                    },
                }
            ]

    with pytest.raises(AttributeError, match="not found"):
        mp.MDVProject._resolve_rows_as_columns_subgroup(
            OneLinkManySubgroups(), "cells", "missing"
        )


def test_resolve_rows_as_columns_no_fallback_when_multiple_single_subgroup_links():
    class TwoLinks:
        def get_links(self, datasource, filter=None):
            del datasource
            return [
                {
                    "datasource": "a",
                    "link": {
                        "rows_as_columns": {
                            "subgroups": {"x": {"name": "a_mat", "type": "dense"}}
                        }
                    },
                },
                {
                    "datasource": "b",
                    "link": {
                        "rows_as_columns": {
                            "subgroups": {"y": {"name": "b_mat", "type": "sparse"}}
                        }
                    },
                },
            ]

    with pytest.raises(AttributeError, match="not found"):
        mp.MDVProject._resolve_rows_as_columns_subgroup(TwoLinks(), "cells", "missing")


def test_get_h5_handle_non_lock_error_reraises_immediately(tmp_path):
    p = mp.MDVProject(str(tmp_path / "proj"), delete_existing=True)
    p.h5file = str(tmp_path / "proj" / "missing.h5")

    with patch.object(
        mp.h5py, "File", side_effect=FileNotFoundError("no file")
    ) as mock_file:
        with pytest.raises(FileNotFoundError):
            p._get_h5_handle(read_only=True)
        assert mock_file.call_count == 1


def test_get_h5_handle_retries_on_lock_contention(tmp_path):
    p = mp.MDVProject(str(tmp_path / "proj"), delete_existing=True)
    p.h5file = str(tmp_path / "proj" / "dat.h5")
    handle = MagicMock()

    with patch.object(
        mp.h5py,
        "File",
        side_effect=[BlockingIOError("unable to lock file"), handle],
    ) as mock_file:
        with patch.object(mp.time, "sleep"):
            out = p._get_h5_handle(read_only=True)
    assert out is handle
    assert mock_file.call_count == 2
