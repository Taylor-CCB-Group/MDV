import json
import threading
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Any

from . import JOBS_DIRNAME
from .jobstore import Status
from .manager import JobManager

logger = logging.getLogger(__name__)

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

    def __init__(
        self,
        # injection seam: any callable returning a manager-shaped object (tests pass fakes)
        manager_factory: Callable[[Any], Any] = JobManager,
        interval: float = 1.0,
    ):
        self._managers: dict[str, JobManager] = {}
        self._manager_factory = manager_factory
        self._wake = threading.Event()
        self._thread: threading.Thread | None = None
        self._interval = interval

    def get_or_create(self, project) -> JobManager:
        if project.id not in self._managers:
            self._managers[project.id] = self._manager_factory(project)
        return self._managers[project.id]

    def tick_all(self) -> None:
        """
        Advance each registered manager once. The driver calls this each cycle
        One manager's failure is logged and skipped so that it doesn't exit the loop for
        all the projects (ADR0012)
        """
        for manager in list(self._managers.values()):
            try:
                manager.tick()
            except Exception:
                logger.exception(
                    "job driver: tick failed for project %s", manager.project.id
                )

    def has_active(self) -> bool:
        """True if any manager holds a non-terminal record"""
        for manager in self._managers.values():
            if any(r.status in _INFLIGHT for r in manager.store.load_all()):
                return True
        return False

    def nudge(self) -> None:
        """Wake the driver at once (submit calls this after writing the queued record)"""
        self._wake.set()

    def _run(self) -> None:
        while True:
            self.tick_all()
            # fast poll while work is in flight so completion is noticed
            # within an interval; block until a nudge while everything is
            # idle (ADR0012)
            if self.has_active():
                self._wake.wait(self._interval)
            else:
                self._wake.wait()
            self._wake.clear()

    def start(self) -> None:
        """Launch the single daemon driver thread, idempotent"""
        if self._thread is None:
            self._thread = threading.Thread(
                target=self._run, name="mdv-job-driver", daemon=True
            )
            self._thread.start()

    def recovery_scan(self, projects) -> list[str]:
        """ADR:0012: at startup, build and reconcile a manager only for projects with in-flight jobs; leave the rest
        for lazy get_or_create. Returns the ids built"""
        built = []
        for project in projects:
            if _has_inflight_records(project):
                self.get_or_create(project)
                built.append(project.id)
        return built
