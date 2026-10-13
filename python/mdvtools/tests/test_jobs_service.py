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
    from typing import Any
    from mdvtools.jobs.service import JobService

    class FakeProject:
        def __init__(self, pid):
            self.id = pid
            self.dir = "/unused"  # tick_all never touches disk with fake managers

    class FakeManager:
        def __init__(self, project, on_submit=None):
            self.project = project
            self.ticks = 0

        def tick(self):
            self.ticks += 1

    service = JobService(manager_factory=FakeManager)
    # the registry is typed to hand back real JobManagers; these are stand-ins
    m1: Any = service.get_or_create(FakeProject("p1"))
    m2: Any = service.get_or_create(FakeProject("p2"))

    service.tick_all()
    assert m1.ticks == 1 and m2.ticks == 1

    service.tick_all()
    assert m1.ticks == 2 and m2.ticks == 2

def test_tick_all_continues_when_one_manager_fails():
    from typing import Any
    from mdvtools.jobs.service import JobService

    class FakeProject:
        def __init__(self, pid):
            self.id = pid
            self.dir = "/unused"

    class BoomManager:
        def __init__(self, project, on_submit=None):
            self.project = project

        def tick(self):
            raise RuntimeError("boom")

    class OkManager:
        def __init__(self, project, on_submit=None):
            self.project = project
            self.ticks = 0

        def tick(self):
            self.ticks += 1

    factories = {"boom": BoomManager, "ok": OkManager}
    service = JobService(
        manager_factory=lambda p, on_submit=None: factories[p.id](p, on_submit=on_submit)
    )
    service.get_or_create(FakeProject("boom"))  # ticks first, raises
    ok: Any = service.get_or_create(FakeProject("ok"))

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

def test_get_or_create_wires_service_nudge_as_on_submit():
    from typing import Any
    from mdvtools.jobs.service import JobService

    class FakeProject:
        def __init__(self, pid):
            self.id = pid
            self.dir = "/unused"

    class FakeManager:
        def __init__(self, project, on_submit=None):
            self.project = project
            self.on_submit = on_submit

    service = JobService(manager_factory=FakeManager)
    mgr: Any = service.get_or_create(FakeProject("p1"))

    # the manager's submit hook is the service's own nudge, so a submit wakes the driver
    assert mgr.on_submit == service.nudge


def test_recovery_scan_survives_a_corrupt_record_and_hands_it_to_quarantine(tmp_path):
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

    # a project whose only record is unparseable
    corrupt = FakeProject("corrupt")
    records = Path(corrupt.dir) / JOBS_DIRNAME / "records"
    records.mkdir(parents=True)
    (records / "broken.json").write_text("{not json")

    # a healthy in-flight project scanned after it
    inflight = FakeProject("inflight")
    JobStore(Path(inflight.dir) / JOBS_DIRNAME).new("concat_columns", {})

    built = JobService().recovery_scan([corrupt, inflight])

    assert built == ["corrupt", "inflight"]  # the scan did not abort
    assert (records / "quarantine" / "broken.json").exists()  # building its manager quarantined it


class _DirProject:
    def __init__(self, tmp_path, pid="proj1"):
        self.id = pid
        d = tmp_path / pid
        d.mkdir()
        self.dir = str(d)


def test_env_without_slurm_offers_only_local(tmp_path):
    from mdvtools.jobs.service import JobService, manager_factory_from_env
    from mdvtools.jobs.executor import LocalSubprocessExecutor
    from mdvtools.jobs.workspace import default_workspace_root

    project = _DirProject(tmp_path)
    service = JobService(manager_factory=manager_factory_from_env({}))
    mgr = service.get_or_create(project)

    assert list(mgr.executors) == ["local"]
    local = mgr.executors["local"]
    assert isinstance(local.executor, LocalSubprocessExecutor)
    assert local.workspace_root == default_workspace_root(project)
    assert local.max_concurrent == 2
    assert mgr.on_submit == service.nudge     # the factory still wires the driver's nudge


def test_env_with_slurm_cli_adds_slurm_after_local(tmp_path):
    from mdvtools.jobs.service import manager_factory_from_env
    from mdvtools.jobs.executor import SlurmExecutor

    project = _DirProject(tmp_path)
    shared = tmp_path / "shared"
    factory = manager_factory_from_env(
        {"MDV_JOBS_SLURM": "cli", "MDV_JOBS_SLURM_WORKSPACE_ROOT": str(shared)}
    )
    mgr = factory(project, on_submit=None)

    assert list(mgr.executors) == ["local", "slurm"]    # local stays the default
    slurm = mgr.executors["slurm"]
    assert isinstance(slurm.executor, SlurmExecutor)
    assert slurm.workspace_root == shared / "proj1"     # one directory per project
    assert slurm.max_concurrent is None                 # Slurm queues the work itself


def test_env_with_slurm_but_no_workspace_root_fails_at_build(tmp_path):
    import pytest
    from mdvtools.jobs.service import manager_factory_from_env

    # raised when the factory is built (server start), not on the first job request
    with pytest.raises(ValueError, match="MDV_JOBS_SLURM_WORKSPACE_ROOT"):
        manager_factory_from_env({"MDV_JOBS_SLURM": "cli"})


def test_env_with_unknown_slurm_mode_fails_at_build(tmp_path):
    import pytest
    from mdvtools.jobs.service import manager_factory_from_env

    with pytest.raises(ValueError, match="sbatch"):
        manager_factory_from_env(
            {"MDV_JOBS_SLURM": "sbatch", "MDV_JOBS_SLURM_WORKSPACE_ROOT": str(tmp_path)}
        )


def test_env_slurm_python_sets_the_worker_python(tmp_path):
    import sys
    from mdvtools.jobs.service import manager_factory_from_env
    from mdvtools.jobs.executor import SlurmExecutor

    settings ={"MDV_JOBS_SLURM": "cli", "MDV_JOBS_SLURM_WORKSPACE_ROOT": str(tmp_path / "shared")}
    default = manager_factory_from_env(settings)(_DirProject(tmp_path, "p1"), on_submit=None)
    chosen = manager_factory_from_env(
        {**settings, "MDV_JOBS_SLURM_PYTHON": "/shared/venv/bin/python"}
    )(_DirProject(tmp_path, "p2"), on_submit=None)

    default_slurm = default.executors["slurm"].executor
    chosen_slurm = chosen.executors["slurm"].executor
    assert isinstance(default_slurm, SlurmExecutor)
    assert isinstance(chosen_slurm, SlurmExecutor)
    # the server's own Python may not exist on the compute nodes
    assert default_slurm.python == sys.executable
    assert chosen_slurm.python == "/shared/venv/bin/python"


def test_env_with_slurm_rest_offers_slurm_through_slurmrestd(tmp_path, fake_slurmrestd):
    from mdvtools.jobs.service import manager_factory_from_env
    from mdvtools.jobs.executor import SlurmRestExecutor

    url, seen = fake_slurmrestd
    token_file = tmp_path / "token"
    token_file.write_text("jwt\n")
    shared = tmp_path / "shared"
    factory = manager_factory_from_env({
        "MDV_JOBS_SLURM": "rest",
        "MDV_JOBS_SLURM_WORKSPACE_ROOT": str(shared),
        "MDV_JOBS_SLURM_PYTHON": "/shared/venv/bin/python",
        "MDV_JOBS_SLURM_REST_URL": url,
        "MDV_JOBS_SLURM_REST_TOKEN_FILE": str(token_file),
        "MDV_JOBS_SLURM_REST_USER": "mdv",
        "MDV_JOBS_SLURM_REST_API_VERSION": "v0.0.42",
    })
    mgr = factory(_DirProject(tmp_path), on_submit=None)

    assert list(mgr.executors) == ["local", "slurm"]    # the user-facing name stays slurm
    slurm = mgr.executors["slurm"]
    assert isinstance(slurm.executor, SlurmRestExecutor)
    assert slurm.workspace_root == shared / "proj1"
    assert slurm.max_concurrent is None

    ws = tmp_path / "ws"
    ws.mkdir()
    slurm.executor.submit("umap_worker:run", ws)

    [req] = seen
    assert req["path"] == "/slurm/v0.0.42/job/submit"
    assert req["token"] == "jwt"
    assert req["user"] == "mdv"
    assert "/shared/venv/bin/python -m mdvtools.jobs.run_worker" in req["body"]["script"]


def test_env_with_slurm_rest_defaults_the_api_version(tmp_path, fake_slurmrestd):
    from mdvtools.jobs.service import manager_factory_from_env

    url, seen = fake_slurmrestd
    token_file = tmp_path / "token"
    token_file.write_text("jwt")
    mgr = manager_factory_from_env({
        "MDV_JOBS_SLURM": "rest",
        "MDV_JOBS_SLURM_WORKSPACE_ROOT": str(tmp_path / "shared"),
        "MDV_JOBS_SLURM_REST_URL": url,
        "MDV_JOBS_SLURM_REST_TOKEN_FILE": str(token_file),
    })(_DirProject(tmp_path), on_submit=None)
    ws = tmp_path / "ws"
    ws.mkdir()

    mgr.executors["slurm"].executor.submit("umap_worker:run", ws)

    assert seen[0]["path"] == "/slurm/v0.0.44/job/submit"
    assert seen[0]["user"] is None      # no user header unless the setting is given


def test_env_with_slurm_rest_but_missing_settings_fails_at_build(tmp_path):
    import pytest
    from mdvtools.jobs.service import manager_factory_from_env

    complete = {
        "MDV_JOBS_SLURM": "rest",
        "MDV_JOBS_SLURM_WORKSPACE_ROOT": str(tmp_path),
        "MDV_JOBS_SLURM_REST_URL": "http://slurmrestd:6820",
        "MDV_JOBS_SLURM_REST_TOKEN_FILE": str(tmp_path / "token"),
    }
    for missing in ("MDV_JOBS_SLURM_REST_URL", "MDV_JOBS_SLURM_REST_TOKEN_FILE"):
        settings = {k: v for k, v in complete.items() if k != missing}
        with pytest.raises(ValueError, match=missing):
            manager_factory_from_env(settings)
