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
