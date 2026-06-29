from mdvtools.jobs.jobstore import JobStore, Status, JobRecord


def test_new_writes_record_to_disk(tmp_path):
    store = JobStore(tmp_path)
    rec = store.new("concat_columns", {"datasource": "cells"})

    assert rec.status == Status.QUEUED.value
    assert (tmp_path / "records" / f"{rec.job_id}.json").exists()
    # reloads from disk as an equal record
    reloaded = {r.job_id: r for r in store.load_all()}
    assert reloaded[rec.job_id].tool_id == "concat_columns"


def test_reconcile_requeues_in_flight_status(tmp_path):
    store = JobStore(tmp_path)
    staging = store.set(store.new("t", {}), Status.STAGING)
    running = store.set(
        store.new("t", {}), Status.RUNNING, handle={"kind": "local", "ref": "123"}
    )
    ingesting = store.set(store.new("t", {}), Status.INGESTING)
    done = store.set(store.new("t", {}), Status.DONE)
    failed = store.set(store.new("t", {}), Status.FAILED)

    store.reconcile_on_boot()

    by_id = {r.job_id: r for r in store.load_all()}
    # all three in-flight status -> queued, handle cleared
    for rec in (staging, running, ingesting):
        assert by_id[rec.job_id].status == Status.QUEUED.value
        assert by_id[rec.job_id].handle is None

    # terminal status remained untouched
    assert by_id[done.job_id].status == Status.DONE.value
    assert by_id[failed.job_id].status == Status.FAILED.value
