import numpy as np
import pytest
from pathlib import Path
from mdvtools.jobs.jobstore import JobStore, Status


def test_new_writes_record_to_disk(tmp_path):
    store = JobStore(tmp_path)
    rec = store.new("concat_columns", {"datasource": "cells"})

    assert rec.status == Status.QUEUED.value
    assert (tmp_path / "records" / f"{rec.job_id}.json").exists()
    # reloads from disk as an equal record
    reloaded = {r.job_id: r for r in store.load_all()}
    assert reloaded[rec.job_id].tool_id == "concat_columns"

def test_set_failed_with_error_persists_message(tmp_path):
    store = JobStore(tmp_path)
    rec = store.new("concat_columns", {"datasource": "cells"})

    store.set(rec, Status.FAILED, error="ingest blew up: column 'x' missing")

    reloaded = {r.job_id: r for r in store.load_all()}[rec.job_id]
    assert reloaded.status == Status.FAILED.value
    assert reloaded.error == "ingest blew up: column 'x' missing"


def test_failed_write_leaves_previous_record_intact(tmp_path, monkeypatch):
    store = JobStore(tmp_path)
    rec = store.new("concat_columns", {"datasource": "cells"})

    real_write_text = Path.write_text

    def write_half_then_crash(self, data, *args, **kwargs):
        real_write_text(self, data[: len(data) // 2], *args, **kwargs)
        raise OSError("disk full")

    monkeypatch.setattr(Path, "write_text", write_half_then_crash)
    with pytest.raises(OSError):
        store.set(rec, Status.RUNNING)
    monkeypatch.undo()

    reloaded = {r.job_id: r for r in store.load_all()}
    assert reloaded[rec.job_id].status == Status.QUEUED.value


@pytest.mark.parametrize(
    "contents",
    ["{not json", '{"job_id": "x", "unknown_field": 1}'],
    ids=["bad_json", "wrong_fields"],
)
def test_load_all_quarantines_malformed_record_and_keeps_the_rest(tmp_path, contents):
    store = JobStore(tmp_path)
    good = store.new("concat_columns", {"datasource": "cells"})
    bad = tmp_path / "records" / "broken.json"
    bad.write_text(contents)

    recs = store.load_all()

    assert [r.job_id for r in recs] == [good.job_id]
    assert not bad.exists()  # out of the active set
    assert (tmp_path / "records" / "quarantine" / "broken.json").read_text() == contents  # kept for inspection


def test_load_all_skips_unreadable_record_and_retries_it_next_pass(tmp_path, monkeypatch):
    store = JobStore(tmp_path)
    good = store.new("concat_columns", {"datasource": "cells"})
    flaky = store.new("concat_columns", {"datasource": "cells"})
    flaky_path = tmp_path / "records" / f"{flaky.job_id}.json"

    real_read_text = Path.read_text

    def locked_once(self, *args, **kwargs):
        if self == flaky_path:
            raise OSError("resource temporarily unavailable")
        return real_read_text(self, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", locked_once)
    first = store.load_all()
    monkeypatch.undo()

    assert [r.job_id for r in first] == [good.job_id]  # skipped this pass
    assert flaky_path.exists()  # transient: not quarantined
    assert {r.job_id for r in store.load_all()} == {good.job_id, flaky.job_id}  # back next pass


def test_new_with_rows_writes_sorted_sidecar_and_records_hash_and_count(tmp_path):
    store = JobStore(tmp_path)
    rec = store.new("umap", {"datasource": "cells"}, rows=[4, 0, 2])

    assert (tmp_path / "selections" / f"{rec.job_id}.npy").exists()
    rows = store.load_rows(rec.job_id)
    assert rows is not None and rows.dtype == np.uint32
    np.testing.assert_array_equal(rows, [0, 2, 4])

    reloaded = {r.job_id: r for r in store.load_all()}[rec.job_id]
    assert reloaded.selection_rows == 3
    assert reloaded.input_filter_hash is not None
    assert reloaded.input_filter_hash == rec.input_filter_hash


def test_same_rows_in_any_order_have_the_same_hash(tmp_path):
    store = JobStore(tmp_path)
    a = store.new("umap", {}, rows=[2, 0])
    b = store.new("umap", {}, rows=[0, 2])
    c = store.new("umap", {}, rows=[0, 3])

    assert a.input_filter_hash == b.input_filter_hash
    assert a.input_filter_hash != c.input_filter_hash


def test_new_without_rows_has_no_selection(tmp_path):
    store = JobStore(tmp_path)
    rec = store.new("umap", {"datasource": "cells"})

    assert rec.selection_rows is None
    assert rec.input_filter_hash is None
    assert store.load_rows(rec.job_id) is None


def test_failed_sidecar_write_leaves_no_record(tmp_path, monkeypatch):
    # the sidecar is written first, so a record never points at rows that are not on disk
    store = JobStore(tmp_path)

    def broken_save(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(np, "save", broken_save)
    with pytest.raises(OSError):
        store.new("umap", {}, rows=[0, 1])

    assert store.load_all() == []
