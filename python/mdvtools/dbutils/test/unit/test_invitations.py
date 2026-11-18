"""The Invitation table that MDV-Admin writes and the sign-in callback reads.

ADR 0019: an Invitation offers one email address access through one Sign-in
method, with a Project role on each chosen project. It lasts 14 days, and a
resend starts a new 14 days. ADR 0021: the email matches ignoring case.
"""

from datetime import datetime

import pytest
from flask import Flask

from mdvtools.dbutils.dbmodels import Invitation, InvitationProject, Project, User, db
from mdvtools.dbutils.dbservice import ProjectService
from mdvtools.dbutils.mdv_server_app import ensure_invitation_tables


@pytest.fixture()
def app():
    app = Flask(__name__)
    app.config.update(
        SQLALCHEMY_DATABASE_URI="sqlite:///:memory:",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
    )
    db.init_app(app)
    with app.app_context():
        db.create_all()
        yield app
        db.session.remove()
        db.drop_all()


@pytest.fixture()
def admin(app):
    user = User(email="admin@example.org", auth_id="auth0|admin", is_active=True, is_admin=True)
    db.session.add(user)
    db.session.commit()
    return user


PASSWORD = "Username-Password-Authentication"
SENT = datetime(2026, 11, 5, 10, 0)


def new_invitation(admin, email="ann@example.org", **fields):
    invitation = Invitation(email=email, connection=PASSWORD, invited_by=admin.id, **fields)
    db.session.add(invitation)
    db.session.commit()
    return invitation


def test_email_is_stored_lowercase_without_surrounding_spaces(admin):
    invitation = new_invitation(admin, email="  Ann.Lee@Example.ORG ")

    assert invitation.email == "ann.lee@example.org"

    invitation.email = "ANN@EXAMPLE.ORG"
    assert invitation.email == "ann@example.org"


def test_invitation_is_open_for_14_days_after_it_is_sent(admin):
    invitation = new_invitation(admin, sent_at=SENT)

    assert invitation.expires_at == datetime(2026, 11, 19, 10, 0)
    assert invitation.is_open(at=datetime(2026, 11, 19, 9, 59))
    assert not invitation.is_open(at=datetime(2026, 11, 19, 10, 0))


def test_resend_starts_a_new_14_days(admin):
    invitation = new_invitation(admin, sent_at=SENT)

    invitation.resend(at=datetime(2026, 11, 25, 9, 0))
    db.session.commit()

    assert invitation.sent_at == datetime(2026, 11, 25, 9, 0)
    assert invitation.expires_at == datetime(2026, 12, 9, 9, 0)
    assert invitation.is_open(at=datetime(2026, 12, 1, 9, 0))


@pytest.mark.parametrize("closed_by", ["accepted_at", "revoked_at"])
def test_accepted_or_revoked_invitation_is_closed(admin, closed_by):
    invitation = new_invitation(admin, sent_at=SENT, **{closed_by: datetime(2026, 11, 6, 10, 0)})

    assert not invitation.is_open(at=datetime(2026, 11, 7, 10, 0))


def test_open_invitations_for_an_email_ignore_case_and_skip_closed_ones(admin):
    now = datetime(2026, 11, 20, 10, 0)
    open_one = new_invitation(admin, sent_at=datetime(2026, 11, 10, 10, 0))
    new_invitation(admin, sent_at=SENT)  # expired on 19 November
    new_invitation(admin, sent_at=datetime(2026, 11, 10, 10, 0), accepted_at=datetime(2026, 11, 11, 10, 0))
    new_invitation(admin, sent_at=datetime(2026, 11, 10, 10, 0), revoked_at=datetime(2026, 11, 11, 10, 0))
    new_invitation(admin, email="bob@example.org", sent_at=datetime(2026, 11, 10, 10, 0))

    assert Invitation.open_for_email(" Ann@Example.org", at=now) == [open_one]


def test_invitation_carries_a_project_role_on_each_chosen_project(admin):
    db.session.add_all([Project(id=1, name="atlas", path="/projects/1"), Project(id=2, name="tumour", path="/projects/2")])
    new_invitation(
        admin,
        projects=[InvitationProject(project_id=1, role="viewer"), InvitationProject(project_id=2, role="owner")],
    )

    stored = Invitation.query.one()
    assert {(p.project_id, p.role) for p in stored.projects} == {(1, "viewer"), (2, "owner")}


def test_project_role_is_viewer_editor_or_owner():
    with pytest.raises(ValueError):
        InvitationProject(project_id=1, role="admin")


def test_purging_a_project_removes_it_from_invitations(admin, tmp_path):
    # Postgres enforces the foreign key, so a leftover row would make the purge fail there.
    (tmp_path / "1").mkdir()
    (tmp_path / "2").mkdir()
    db.session.add_all([
        Project(id=1, name="atlas", path=str(tmp_path / "1"), is_deleted=True),
        Project(id=2, name="tumour", path=str(tmp_path / "2")),
    ])
    new_invitation(
        admin,
        projects=[InvitationProject(project_id=1, role="viewer"), InvitationProject(project_id=2, role="editor")],
    )

    ok, error = ProjectService.purge_deleted_project(1)

    assert ok is True, error
    assert [(p.project_id, p.role) for p in Invitation.query.one().projects] == [(2, "editor")]


def test_startup_adds_the_invitation_tables_to_an_existing_database(tmp_path):
    # create_all() runs only on an empty database, so an established deployment
    # gets these tables from this step alone.
    app = Flask(__name__)
    app.config.update(
        SQLALCHEMY_DATABASE_URI=f"sqlite:///{tmp_path / 'existing.sqlite3'}",
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
    )
    db.init_app(app)
    with app.app_context():
        db.metadata.create_all(db.engine, tables=[User.__table__, Project.__table__])

        ensure_invitation_tables()

        tables = set(db.inspect(db.engine).get_table_names())
        db.session.remove()
    assert {"invitations", "invitation_projects"} <= tables
