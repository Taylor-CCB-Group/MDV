"""Tests for usage telemetry recording and its interaction with project deletion."""

import random
from datetime import datetime

import pytest
from flask import Flask

from mdvtools.dbutils.dbmodels import Project, UsageEvent, User, db
from mdvtools.dbutils import dbservice
from mdvtools.dbutils.dbservice import ProjectService, UsageEventService
from mdvtools.dbutils.seed_usage_events import clear_seeded, generate, summarise


@pytest.fixture()
def app(tmp_path):
    app = Flask(__name__)
    app.config.update(
        SECRET_KEY="test-only",
        SQLALCHEMY_DATABASE_URI="sqlite:///:memory:",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
        ENABLE_AUTH=False,
        projects_base_dir=str(tmp_path),
    )
    db.init_app(app)
    with app.app_context():
        db.create_all()
        yield app
        db.session.remove()
        db.drop_all()


def add_user(user_id=7):
    user = User(id=user_id, email=f"user{user_id}@example.org", auth_id=f"auth0|{user_id}")
    db.session.add(user)
    db.session.commit()
    return user


def add_project(tmp_path, project_id=3, *, deleted=False):
    path = tmp_path / str(project_id)
    path.mkdir()
    project = Project(
        id=project_id,
        name=f"project-{project_id}",
        path=str(path),
        is_deleted=deleted,
        access_level="editable",
    )
    db.session.add(project)
    db.session.commit()
    return project


class TestRecording:
    def test_login_records_one_row_with_no_project(self, app):
        add_user()
        assert UsageEventService.record_login(7) is True

        events = UsageEvent.query.all()
        assert len(events) == 1
        assert events[0].event_type == "login"
        assert events[0].project_id is None
        assert events[0].view_name is None
        assert events[0].occurred_at is not None

    def test_project_open_records_the_project(self, app, tmp_path):
        add_user()
        add_project(tmp_path)
        UsageEventService.record_project_open(7, 3)

        event = UsageEvent.query.first()
        assert event.event_type == "project_open"
        assert event.project_id == 3
        assert event.view_name is None

    def test_view_open_records_the_view_name(self, app, tmp_path):
        add_user()
        add_project(tmp_path)
        UsageEventService.record_view_open(7, 3, "QC plots")

        event = UsageEvent.query.first()
        assert event.event_type == "view_open"
        assert event.project_id == 3
        assert event.view_name == "QC plots"

    def test_an_overlong_view_name_is_truncated_not_rejected(self, app, tmp_path):
        """View names arrive from an untrusted request body."""
        add_user()
        add_project(tmp_path)
        assert UsageEventService.record_view_open(7, 3, "x" * 500) is True
        assert len(UsageEvent.query.first().view_name) == 128

    def test_a_missing_user_id_records_nothing_and_does_not_raise(self, app):
        assert UsageEventService.record_event(None, "login") is False
        assert UsageEvent.query.count() == 0

    def test_a_failure_is_swallowed_rather_than_raised(self, app, monkeypatch):
        """Telemetry must never turn a working page into an error."""
        def explode(*args, **kwargs):
            raise RuntimeError("database is on fire")

        monkeypatch.setattr(db.engine, "begin", explode)
        assert UsageEventService.record_login(7) is False

    def test_recording_does_not_disturb_the_caller_uncommitted_work(self, app, tmp_path):
        """Writes go via db.engine.begin(), not db.session, on purpose.

        The login hook fires straight after a user-activation commit, so rolling
        back the shared session here would discard the caller's work.
        """
        add_user()
        pending = Project(id=99, name="pending", path=str(tmp_path / "pending"))
        db.session.add(pending)

        UsageEventService.record_login(7)

        db.session.commit()
        assert Project.query.get(99) is not None


class TestCaptureGuards:
    """record_usage_event is the guard layer between MDV routes and the service."""

    @staticmethod
    def _call(app, project_id, *, backend_db=True, enable_auth=True, user=None, view=None):
        from flask import session as flask_session

        from mdvtools.server import record_usage_event

        app.config["ENABLE_AUTH"] = enable_auth
        with app.test_request_context("/"):
            if user is not None:
                flask_session["user"] = user
            return record_usage_event(project_id, "project_open", backend_db, view_name=view)

    def test_a_string_project_id_is_coerced_to_the_row_id(self, app, tmp_path):
        """MDVProject.id is built as str(project.id), so this is the normal case."""
        add_user()
        add_project(tmp_path)

        assert self._call(app, "3", user={"id": 7}) is True

        event = UsageEvent.query.first()
        assert event.project_id == 3

    def test_a_non_numeric_project_id_records_nothing(self, app):
        """Single-project mode uses a directory name, which is not a row id."""
        add_user()
        assert self._call(app, "my-project-folder", user={"id": 7}) is False
        assert UsageEvent.query.count() == 0

    def test_nothing_is_recorded_without_a_backend_database(self, app, tmp_path):
        add_user()
        add_project(tmp_path)
        assert self._call(app, "3", backend_db=False, user={"id": 7}) is False
        assert UsageEvent.query.count() == 0

    def test_nothing_is_recorded_when_auth_is_disabled(self, app, tmp_path):
        add_user()
        add_project(tmp_path)
        assert self._call(app, "3", enable_auth=False, user={"id": 7}) is False
        assert UsageEvent.query.count() == 0

    def test_nothing_is_recorded_without_a_session_user(self, app, tmp_path):
        add_user()
        add_project(tmp_path)
        assert self._call(app, "3", user=None) is False
        assert UsageEvent.query.count() == 0


class TestViewCreation:
    """A view created in the browser is never fetched, so /get_view never fires."""

    def test_creating_a_view_is_recorded_as_its_own_event_type(self, app, tmp_path):
        add_user()
        add_project(tmp_path)
        UsageEventService.record_event(7, "view_create", project_id=3, view_name="New view")

        event = UsageEvent.query.first()
        assert event.event_type == "view_create"
        assert event.view_name == "New view"

    def test_creation_does_not_inflate_the_open_count(self, app, tmp_path):
        """Counting a creation as an open would make every view read one too many."""
        add_user()
        add_project(tmp_path)
        UsageEventService.record_event(7, "view_create", project_id=3, view_name="New view")
        UsageEventService.record_view_open(7, 3, "New view")

        opens = UsageEvent.query.filter_by(event_type="view_open").count()
        creates = UsageEvent.query.filter_by(event_type="view_create").count()
        assert opens == 1
        assert creates == 1


class TestOffSwitch:
    """ENABLE_USAGE_TRACKING lets a deployment record nothing at all."""

    @pytest.fixture()
    def tracking_off(self, monkeypatch):
        monkeypatch.setattr(dbservice, "ENABLE_USAGE_TRACKING", False)

    def test_recording_is_on_by_default(self, app):
        """A deployment that says nothing gets usage recording."""
        assert dbservice.ENABLE_USAGE_TRACKING is True

    def test_nothing_is_recorded_when_switched_off(self, app, tracking_off):
        add_user()
        assert UsageEventService.record_login(7) is False
        assert UsageEvent.query.count() == 0

    def test_the_switch_covers_every_event_type(self, app, tmp_path, tracking_off):
        """Checked in record_event, so a new event type cannot bypass it."""
        add_user()
        add_project(tmp_path)
        UsageEventService.record_login(7)
        UsageEventService.record_project_open(7, 3)
        UsageEventService.record_view_open(7, 3, "Overview")
        UsageEventService.record_event(7, "view_create", project_id=3, view_name="New")
        assert UsageEvent.query.count() == 0

    def test_switching_it_off_does_not_delete_what_was_already_recorded(self, app, tmp_path, monkeypatch):
        """It stops collection. Removing history is retention's job, not this."""
        add_user()
        add_project(tmp_path)
        UsageEventService.record_project_open(7, 3)
        monkeypatch.setattr(dbservice, "ENABLE_USAGE_TRACKING", False)
        UsageEventService.record_project_open(7, 3)
        assert UsageEvent.query.count() == 1


class TestProjectDeletion:
    def test_soft_delete_leaves_usage_history_untouched(self, app, tmp_path):
        add_user()
        add_project(tmp_path)
        UsageEventService.record_project_open(7, 3)

        ProjectService.soft_delete_project(3)

        assert UsageEvent.query.count() == 1
        assert UsageEvent.query.first().project_id == 3

    def test_restoring_a_project_keeps_its_history_attached(self, app, tmp_path):
        add_user()
        add_project(tmp_path)
        UsageEventService.record_project_open(7, 3)

        ProjectService.soft_delete_project(3)
        ProjectService.restore_deleted_project(3)

        assert UsageEvent.query.first().project_id == 3

    def test_purge_succeeds_and_detaches_history_instead_of_deleting_it(self, app, tmp_path):
        """Without this, the foreign key would make emptying the recycle bin fail."""
        add_user()
        add_project(tmp_path, deleted=True)
        UsageEventService.record_view_open(7, 3, "Overview")

        ok, error = ProjectService.purge_deleted_project(3)

        assert ok is True, error
        assert Project.query.get(3) is None

        event = UsageEvent.query.first()
        assert event is not None, "usage history should survive a purge"
        assert event.project_id is None
        assert event.view_name == "Overview"
        assert event.details == {"project_name": "project-3"}


class TestSeedGeneration:
    """generate() is pure, so the shape of demo data can be checked without a database."""

    @staticmethod
    def _users(count=5):
        return [User(id=i, email=f"u{i}@example.org", auth_id=f"a{i}") for i in range(1, count + 1)]

    @staticmethod
    def _projects(count=3):
        return [Project(id=i, name=f"p{i}", path=f"/tmp/p{i}") for i in range(1, count + 1)]

    def test_it_produces_events(self):
        events = generate(self._users(), self._projects(), 30, random.Random(1))
        assert len(events) > 0

    def test_a_session_always_starts_with_a_sign_in(self):
        events = generate(self._users(), self._projects(), 30, random.Random(1))
        assert events[0][2] == "login"

    def test_a_sign_in_has_no_project(self):
        events = generate(self._users(), self._projects(), 30, random.Random(1))
        assert all(e[1] is None for e in events if e[2] == "login")

    def test_view_events_always_name_a_view(self):
        events = generate(self._users(), self._projects(), 30, random.Random(1))
        assert all(e[3] for e in events if e[2] in ("view_open", "view_create"))

    def test_project_opens_never_name_a_view(self):
        events = generate(self._users(), self._projects(), 30, random.Random(1))
        assert all(e[3] is None for e in events if e[2] == "project_open")

    def test_nothing_is_dated_in_the_future(self):
        events = generate(self._users(), self._projects(), 30, random.Random(1))
        assert all(e[4] <= datetime.now() for e in events)

    def test_the_same_seed_gives_the_same_data(self):
        a = generate(self._users(), self._projects(), 30, random.Random(7))
        b = generate(self._users(), self._projects(), 30, random.Random(7))
        assert a == b

    def test_weekends_are_quieter_than_weekdays(self):
        events = generate(self._users(), self._projects(), 60, random.Random(3))
        weekend = sum(1 for e in events if e[4].weekday() >= 5)
        weekday = sum(1 for e in events if e[4].weekday() < 5)
        assert weekend < weekday / 2, "weekend dip is what makes the chart read as real"

    def test_activity_is_uneven_across_people(self):
        """A demo needs the spread - somebody heavy, somebody who never shows up."""
        events = generate(self._users(), self._projects(), 60, random.Random(5))
        per_user = {}
        for user_id, *_ in events:
            per_user[user_id] = per_user.get(user_id, 0) + 1
        assert len(set(per_user.values())) > 1
        assert len(per_user) < 5, "at least one seeded user should be dormant"


class TestSeedMarker:
    """--clear must remove only what the script wrote, however often it has run."""

    def test_clear_removes_seeded_rows_and_leaves_real_ones(self, app, tmp_path):
        add_user()
        add_project(tmp_path)
        UsageEventService.record_project_open(7, 3)
        db.session.add(UsageEvent(
            user_id=7, project_id=3, event_type="view_open",
            view_name="seeded", occurred_at=datetime.now(), details={"seeded": True},
        ))
        db.session.commit()
        assert UsageEvent.query.count() == 2

        assert clear_seeded() == 1

        remaining = UsageEvent.query.all()
        assert len(remaining) == 1
        assert remaining[0].event_type == "project_open"

    def test_clear_is_safe_to_run_when_nothing_was_seeded(self, app, tmp_path):
        add_user()
        add_project(tmp_path)
        UsageEventService.record_project_open(7, 3)
        assert clear_seeded() == 0
        assert UsageEvent.query.count() == 1

    def test_clearing_twice_is_harmless(self, app, tmp_path):
        """--reset calls this before seeding, so it runs on already-clean tables."""
        add_user()
        add_project(tmp_path)
        db.session.add(UsageEvent(
            user_id=7, project_id=3, event_type="login",
            occurred_at=datetime.now(), details={"seeded": True},
        ))
        db.session.commit()
        assert clear_seeded() == 1
        assert clear_seeded() == 0

    def test_a_row_with_other_details_is_not_treated_as_seeded(self, app, tmp_path):
        """A purged project writes project_name into details - it must survive."""
        add_user()
        add_project(tmp_path)
        db.session.add(UsageEvent(
            user_id=7, project_id=None, event_type="project_open",
            occurred_at=datetime.now(), details={"project_name": "Old pilot"},
        ))
        db.session.commit()
        assert clear_seeded() == 0
        assert UsageEvent.query.count() == 1

    def test_summary_separates_seeded_from_real(self, app, tmp_path):
        add_user()
        add_project(tmp_path)
        UsageEventService.record_project_open(7, 3)
        db.session.add(UsageEvent(
            user_id=7, project_id=3, event_type="login",
            occurred_at=datetime.now(), details={"seeded": True},
        ))
        db.session.commit()
        assert summarise() == (2, 1, 1)
