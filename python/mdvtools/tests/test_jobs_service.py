def test_job_service_builds_one_manager_per_project(tmp_path):
    from mdvtools.jobs.service import JobService
    from mdvtools.jobs.manager import JobManager

    class FakeProject:
        def __init__(self, pid):
            self.id = pid
            d = tmp_path / pid
            d.mkdir()
            self.dir = str(d)

    p1 = FakeProject("proj1")
    p2 = FakeProject("proj2")
    service = JobService()

    m1 = service.get_or_create(p1)
    assert isinstance(m1, JobManager)
    assert m1.project is p1

    # same project id returns the same cached manager, not a rebuilt one
    assert service.get_or_create(p1) is m1

    # a different project gets its own manager
    m2 = service.get_or_create(p2)
    assert m2 is not m1
    assert m2.project is p2

def test_recovery_scan_builds_managers_only_for_inflight_projects(tmp_path):
    from pathlib import Path
    from mdvtools.jobs.service import JobService
    from mdvtools.jobs.jobstore import JobStore, Status
    from mdvtools.jobs import JOBS_DIRNAME

    class FakeProject:
        def __init__(self, pid):
            self.id = pid
            d = tmp_path / pid
            d.mkdir()
            self.dir = str(d)

    # a QUEUED record that outlived a restart -> still needs driving
    inflight = FakeProject("inflight")
    JobStore(Path(inflight.dir) / JOBS_DIRNAME).new("concat_columns", {})

    # no jobs/ dir at all -> the cheap stat skips it
    idle = FakeProject("idle")

    # only a terminal record -> nothing to recover
    done = FakeProject("done")
    store = JobStore(Path(done.dir) / JOBS_DIRNAME)
    store.set(store.new("concat_columns", {}), Status.DONE)

    service = JobService()
    built = service.recovery_scan([inflight, idle, done])

    assert built == ["inflight"]
