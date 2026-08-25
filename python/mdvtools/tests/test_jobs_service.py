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

def test_tick_all_ticks_every_manager_once():
    from mdvtools.jobs.service import JobService

    class FakeProject:
        def __init__(self, pid):
            self.id = pid
            self.dir = "/unused"  # tick_all never touches disk with fake managers

    class FakeManager:
        def __init__(self, project):
            self.project = project
            self.ticks = 0

        def tick(self):
            self.ticks += 1

    service = JobService(manager_factory=FakeManager)
    m1 = service.get_or_create(FakeProject("p1"))
    m2 = service.get_or_create(FakeProject("p2"))

    service.tick_all()
    assert m1.ticks == 1 and m2.ticks == 1

    service.tick_all()
    assert m1.ticks == 2 and m2.ticks == 2

def test_tick_all_continues_when_one_manager_fails():
    from mdvtools.jobs.service import JobService

    class FakeProject:
        def __init__(self, pid):
            self.id = pid
            self.dir = "/unused"

    class BoomManager:
        def __init__(self, project):
            self.project = project

        def tick(self):
            raise RuntimeError("boom")

    class OkManager:
        def __init__(self, project):
            self.project = project
            self.ticks = 0

        def tick(self):
            self.ticks += 1

    factories = {"boom": BoomManager, "ok": OkManager}
    service = JobService(manager_factory=lambda p: factories[p.id](p))
    service.get_or_create(FakeProject("boom"))  # ticks first, raises
    ok = service.get_or_create(FakeProject("ok"))

    service.tick_all()  # must not raise

    assert ok.ticks == 1

def test_has_active_true_only_with_inflight_records(tmp_path):
    from pathlib import Path
    from mdvtools.jobs.service import JobService
    from mdvtools.jobs.jobstore import JobStore
    from mdvtools.jobs import JOBS_DIRNAME

    class FakeProject:
        def __init__(self, pid):
            self.id = pid
            d = tmp_path / pid
            d.mkdir()
            self.dir = str(d)

    service = JobService()

    # a project with no records -> nothing in flight
    service.get_or_create(FakeProject("empty"))
    assert service.has_active() is False

    # a project with a QUEUED record -> in flight
    active = FakeProject("active")
    JobStore(Path(active.dir) / JOBS_DIRNAME).new("concat_columns", {})
    service.get_or_create(active)
    assert service.has_active() is True
