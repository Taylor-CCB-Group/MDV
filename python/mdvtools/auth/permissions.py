"""Project roles, read from the database each time they are checked.

Every check queries the user_projects table, so a role written by any process
applies on the next request: this app, MDV-Admin on a shared database, or a
script run with `docker exec`. No copy is kept between requests, so there is
nothing to rebuild or keep in step with the database.
"""

from typing import Iterable, Optional, TypedDict

from mdvtools.dbutils.dbmodels import User, UserProject, db


class ProjectRole(TypedDict):
    can_read: bool
    can_write: bool
    is_owner: bool


def _role(row: UserProject) -> ProjectRole:
    return {
        "can_read": bool(row.can_read),
        "can_write": bool(row.can_write),
        "is_owner": bool(row.is_owner),
    }


def project_role(user_id: int, project_id: int) -> Optional[ProjectRole]:
    """The user's role on one project, or None when they hold no role on it."""
    row = UserProject.query.filter_by(user_id=user_id, project_id=project_id).first()
    return _role(row) if row else None


def project_roles(user_id: int) -> dict[int, ProjectRole]:
    """The user's role on every project they hold one on, keyed by project ID."""
    rows = UserProject.query.filter_by(user_id=user_id).all()
    return {row.project_id: _role(row) for row in rows}


def can_view(role: Optional[ProjectRole]) -> bool:
    """Whether a role lets its holder open the project."""
    return bool(role) and bool(role["can_read"] or role["can_write"] or role["is_owner"])


def project_owner_emails(project_ids: Iterable[int]) -> dict[int, list[str]]:
    """The sorted email addresses of each project's owners, keyed by project ID.

    A project with no owner has no key.
    """
    ids = list(project_ids)
    if not ids:
        return {}
    rows = (
        db.session.query(UserProject.project_id, User.email)
        .join(User, User.id == UserProject.user_id)
        .filter(UserProject.is_owner, UserProject.project_id.in_(ids))
        .all()
    )
    owners: dict[int, list[str]] = {}
    for project_id, email in rows:
        owners.setdefault(project_id, []).append(email)
    for emails in owners.values():
        emails.sort()
    return owners
