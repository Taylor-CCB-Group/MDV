import os
from typing import Optional

from celery import Celery


def make_celery(app=None) -> Celery:
    """Create a Celery app configured from Flask app config or environment.

    Looks for:
      - CELERY_BROKER_URL
      - CELERY_RESULT_BACKEND
    with sensible defaults if not provided.
    """
    broker_url = None
    result_backend = None

    if app is not None:
        broker_url = app.config.get("CELERY_BROKER_URL")
        result_backend = app.config.get("CELERY_RESULT_BACKEND")

    broker_url = broker_url or os.getenv("CELERY_BROKER_URL") or "redis://localhost:6379/0"
    result_backend = result_backend or os.getenv("CELERY_RESULT_BACKEND") or broker_url

    celery = Celery(
        app.import_name if app else __name__,
        broker=broker_url,
        backend=result_backend,
        include=[
            "mdvtools.tasks.sample_tasks",
        ],
    )

    # Namespace CELERY_ for configuration keys to mirror Flask patterns
    celery.conf.update(
        task_track_started=True,
        task_serializer="json",
        result_serializer="json",
        accept_content=["json"],
        timezone=os.getenv("CELERY_TIMEZONE", "UTC"),
        enable_utc=True,
    )

    return celery


