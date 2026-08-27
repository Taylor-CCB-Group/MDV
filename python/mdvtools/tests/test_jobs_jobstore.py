from mdvtools.jobs.jobstore import JobStore, Status


def test_new_writes_record_to_disk(tmp_path):
    store = JobStore(tmp_path)
    rec = store.new("concat_columns", {"datasource": "cells"})

    assert rec.status == Status.QUEUED.value
    assert (tmp_path / "records" / f"{rec.job_id}.json").exists()
    # reloads from disk as an equal record
    reloaded = {r.job_id: r for r in store.load_all()}
    assert reloaded[rec.job_id].tool_id == "concat_columns"
