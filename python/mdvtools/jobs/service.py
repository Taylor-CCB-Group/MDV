from .manager import JobManager

class JobService:
    """
        Process-wide registry of one JobManager per project

        Managers are built lazily on first request and cached by project-id, so the
        server process holds exactly one owner side manager per project
    """

    def __init__(self):
        self._managers: dict[str, JobManager] = {}

    def get_or_create(self, project) -> JobManager:
        if project.id not in self._managers:
            self._managers[project.id] = JobManager(project)
        return self._managers[project.id]
