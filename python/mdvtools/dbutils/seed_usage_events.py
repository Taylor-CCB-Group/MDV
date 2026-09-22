"""Seed realistic usage events so the Admin usage page can be demonstrated.

**Demo and development only.** This writes invented activity into usage_events.
Never run it against a deployment whose numbers anyone relies on.

Every row it writes carries ``details = {"seeded": true}``. That marker is what
makes ``--clear`` safe: it deletes exactly the rows this script created and can
never remove a real event, however many times it has been run.

Examples::

    python -m mdvtools.dbutils.seed_usage_events --days 60 --yes
    python -m mdvtools.dbutils.seed_usage_events --reset --yes
    python -m mdvtools.dbutils.seed_usage_events --clear --yes
    python -m mdvtools.dbutils.seed_usage_events --summary
"""

import argparse
import random
import sys
from datetime import datetime, timedelta

from mdvtools.dbutils.dbmodels import Project, UsageEvent, User, db
from mdvtools.logging_config import get_logger

logger = get_logger(__name__)

SEED_MARKER = {"seeded": True}

# Plausible view names. Real projects have their own, but a demo only needs
# names that read like something a researcher would make.
VIEW_NAMES = [
    "default",
    "Overview",
    "QC plots",
    "Cell types",
    "UMAP",
    "Marker genes",
    "Cluster comparison",
    "Spatial view",
]

# How often someone shows up, as sessions per working day. A demo needs the
# spread: somebody living in the app, somebody who signed up and vanished.
ACTIVITY_LEVELS = [
    ("heavy", 1.6),
    ("regular", 0.7),
    ("occasional", 0.25),
    ("rare", 0.08),
    ("dormant", 0.0),
]


def _session_events(user_id, project, when, rng):
    """One visit: sign in, open a project, look at some views, occasionally make one.

    Mirrors what MDV actually records - a login once per session, a project_open
    for the page, then a view_open per view loaded - so the seeded data exercises
    the same shapes as real data rather than an idealised version of them.
    """
    events = [(user_id, None, "login", None, when)]
    moment = when + timedelta(seconds=rng.randint(5, 40))
    events.append((user_id, project.id, "project_open", None, moment))

    views = rng.sample(VIEW_NAMES, k=min(rng.randint(1, 4), len(VIEW_NAMES)))
    for view in views:
        moment += timedelta(seconds=rng.randint(20, 900))
        events.append((user_id, project.id, "view_open", view, moment))

    if rng.random() < 0.08:
        moment += timedelta(seconds=rng.randint(60, 600))
        events.append((user_id, project.id, "view_create", f"{rng.choice(VIEW_NAMES)} v2", moment))

    return events


def generate(users, projects, days, rng):
    """Build the event tuples without touching the database.

    Separated from the insert so the shape can be checked - and tested - without
    a database, and so a dry run is genuinely free of side effects.
    """
    levels = {user.id: ACTIVITY_LEVELS[index % len(ACTIVITY_LEVELS)] for index, user in enumerate(users)}
    # Give the newest project a head start so "most used" is not a dead heat.
    weights = [max(1, len(projects) - index) for index in range(len(projects))]

    end = datetime.now().replace(hour=0, minute=0, second=0, microsecond=0)
    events = []
    for day_offset in range(days, -1, -1):
        day = end - timedelta(days=day_offset)
        weekend = day.weekday() >= 5
        for user in users:
            _, rate = levels[user.id]
            if weekend:
                rate *= 0.2
            sessions = int(rate) + (1 if rng.random() < (rate % 1) else 0)
            for _ in range(sessions):
                start = day + timedelta(
                    hours=rng.randint(8, 18), minutes=rng.randint(0, 59), seconds=rng.randint(0, 59)
                )
                if start > datetime.now():
                    continue
                project = rng.choices(projects, weights=weights, k=1)[0]
                events.extend(_session_events(user.id, project, start, rng))
    return events


def clear_seeded():
    """Remove only rows this script wrote, identified by the details marker."""
    rows = UsageEvent.query.filter(
        UsageEvent.details.isnot(None), UsageEvent.event_type.isnot(None)
    ).all()
    removed = 0
    for row in rows:
        if isinstance(row.details, dict) and row.details.get("seeded") is True:
            db.session.delete(row)
            removed += 1
    db.session.commit()
    return removed


def summarise():
    total = UsageEvent.query.count()
    seeded = sum(
        1
        for row in UsageEvent.query.all()
        if isinstance(row.details, dict) and row.details.get("seeded") is True
    )
    return total, seeded, total - seeded


def main():
    parser = argparse.ArgumentParser(
        description="Seed demo usage events. Demo and development use only.",
    )
    parser.add_argument("--days", type=int, default=60, help="how far back to generate (default 60)")
    parser.add_argument("--seed", type=int, default=None, help="random seed, for repeatable data")
    parser.add_argument("--clear", action="store_true", help="delete previously seeded rows and stop")
    parser.add_argument(
        "--reset",
        action="store_true",
        help="clear previously seeded rows first, so re-running replaces rather than doubles",
    )
    parser.add_argument("--summary", action="store_true", help="report what is in the table and stop")
    parser.add_argument("--dry-run", action="store_true", help="report what would be written")
    parser.add_argument(
        "--yes",
        action="store_true",
        help="required to write anything - this is fabricated data",
    )
    args = parser.parse_args()

    from mdvtools.dbutils.mdv_server_app import app

    with app.app_context():
        if args.summary:
            total, seeded, real = summarise()
            print(f"usage_events: {total} rows - {seeded} seeded, {real} real")
            return 0

        if args.clear:
            if not args.yes:
                print("Refusing to delete without --yes.")
                return 1
            print(f"Removed {clear_seeded()} seeded row(s). Real events untouched.")
            return 0

        users = User.query.all()
        projects = Project.query.all()
        if not users or not projects:
            print("Nothing to seed against: this database has no users or no projects.")
            return 1

        rng = random.Random(args.seed)
        events = generate(users, projects, args.days, rng)

        if args.dry_run:
            print(f"Would write {len(events)} event(s) across {len(users)} user(s), "
                  f"{len(projects)} project(s), {args.days} day(s).")
            return 0

        if not args.yes:
            print("This writes fabricated activity into usage_events.")
            print("Never run it where the numbers matter. Re-run with --yes to proceed.")
            return 1

        if args.reset:
            print(f"Removed {clear_seeded()} previously seeded row(s).")

        for user_id, project_id, event_type, view_name, occurred_at in events:
            db.session.add(
                UsageEvent(
                    user_id=user_id,
                    project_id=project_id,
                    event_type=event_type,
                    view_name=view_name,
                    occurred_at=occurred_at,
                    details=dict(SEED_MARKER),
                )
            )
        db.session.commit()
        print(f"Wrote {len(events)} seeded event(s) across {args.days} day(s).")
        print("Remove them again with: --clear --yes")
        return 0


if __name__ == "__main__":
    sys.exit(main())
