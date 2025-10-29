from time import sleep
from typing import Dict

from celery import shared_task


@shared_task(name="mdv.sample.add")
def add(x: int, y: int) -> int:
    return x + y


@shared_task(name="mdv.sample.long_running")
def long_running(params: Dict[str, str]) -> Dict[str, str]:
    """Example long-running task to demo progress/state."""
    total_steps = 5
    for i in range(total_steps):
        sleep(1)
    return {"status": "done", "received": params}


