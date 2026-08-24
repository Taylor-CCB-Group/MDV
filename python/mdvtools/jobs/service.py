import json
from pathlib import Path

from . import JOBS_DIRNAME
from .jobstore import Status
from .manager import JobManager

# non-terminal status: a record in any of these still needs the driver to advance it.
_INFLIGHT = frozenset({
    Status.QUEUED.value, # queued is added so that a write-ahead intent isn't lost when the manager restarts
    Status.STAGING.value,
    Status.RUNNING.value,
    Status.INGESTING.value,
})

def _has_inflight_records(project) -> bool:
    """Peek at a project's records without building a manager of its JobStore"""
    records_dir = Path(project.dir) / JOBS_DIRNAME / "records"
    if not records_dir.exists():
        return False
    for p in records_dir.glob("*.json"):
        if json.loads(p.read_text()).get("status") in _INFLIGHT:
            return True
    return False

class JobService:
    """
        Process-wide registry of one JobManager per project

        Managers are built lazily on first request and cached by project-id, so the
        server process holds exactly one owner side manager per project
    """

    def __init__(self, manager_factory=JobManager):
        self._managers: dict[str, JobManager] = {}
        self._manager_factory = manager_factory

    def get_or_create(self, project) -> JobManager:
        if project.id not in self._managers:
            self._managers[project.id] = self._manager_factory(project)
        return self._managers[project.id]

    def tick_all(self) -> None:
        """Advance each registered manager once. The driver calls this each cycle"""
        for manager in list(self._managers.values()):
            manager.tick()

    def recovery_scan(self, projects) -> list[str]:
        """ADR:0012: at startup, build and reconcile a manager only for projects with in-flight jobs; leave the rest
        for lazy get_or_create. Returns the ids built"""
        built = []
        for project in projects:
            if _has_inflight_records(project):
                self.get_or_create(project)
                built.append(project.id)
        return built
