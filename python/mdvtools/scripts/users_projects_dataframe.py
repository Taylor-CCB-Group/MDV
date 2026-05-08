#!/usr/bin/env python3
"""
Collect users and projects from the database and merge them into a DataFrame.

Each row represents one user–project association (from user_projects), with
user and project attributes joined. Output can be printed, saved to CSV, or
returned for use in other code.

Requires database credentials (env or Docker secrets) as for the main MDV app.
"""
import argparse
import sys
from pathlib import Path
from typing import Dict, Any

# Add the parent directory to sys.path to import mdvtools modules
_scripts_dir = Path(__file__).resolve().parent
_root = _scripts_dir.parent.parent
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

import pandas as pd
from auth0.authentication import GetToken  # pyright: ignore[reportMissingModuleSource]
from auth0.management import Auth0  # pyright: ignore[reportMissingModuleSource]

from mdvtools.dbutils.dbmodels import db, User, Project, UserProject
from mdvtools.dbutils.safe_mdv_app import app


def _fetch_auth0_login_stats() -> Dict[str, Dict[str, Any]]:
    """
    Fetch Auth0 login stats keyed by Auth0 user id.

    Returns:
        Dict[auth0_user_id, {"last_login": <iso str|None>, "logins_count": <int|None>}]
    """
    auth0_domain = app.config.get("AUTH0_DOMAIN")
    client_id = app.config.get("AUTH0_CLIENT_ID")
    client_secret = app.config.get("AUTH0_CLIENT_SECRET")
    auth0_db_connection = app.config.get("AUTH0_DB_CONNECTION")

    if not all([auth0_domain, client_id, client_secret, auth0_db_connection]):
        raise RuntimeError(
            "Missing Auth0 configuration. Expected AUTH0_DOMAIN, AUTH0_CLIENT_ID, "
            "AUTH0_CLIENT_SECRET, and AUTH0_DB_CONNECTION in app config."
        )

    audience = f"https://{auth0_domain}/api/v2/"
    token_client = GetToken(
        domain=auth0_domain,
        client_id=client_id,
        client_secret=client_secret,
    )
    mgmt_api_token = token_client.client_credentials(audience=audience)["access_token"]
    auth0 = Auth0(auth0_domain, mgmt_api_token)

    result: Dict[str, Dict[str, Any]] = {}
    page = 0
    per_page = 100
    while True:
        users_page = auth0.users.list(
            q=f'identities.connection:"{auth0_db_connection}"',
            page=page,
            per_page=per_page,
        )
        user_list = users_page.get("users", [])
        if not user_list:
            break

        for user in user_list:
            auth0_id = user.get("user_id")
            if not auth0_id:
                continue
            result[auth0_id] = {
                "last_login": user.get("last_login"),
                "logins_count": user.get("logins_count"),
            }

        if len(user_list) < per_page:
            break
        page += 1

    return result


def load_users_projects_df(include_auth0_login_stats: bool = False):
    """
    Load users and projects from the database and merge into a single DataFrame.

    Returns a DataFrame with one row per user–project association, with columns:
    - user_id, email, first_name, last_name, auth_id, is_active, is_admin
    - project_id, project_name, project_path, project_owner_id, is_public, is_deleted
    - can_read, can_write, is_owner
    """
    with app.app_context():
        # Query all user–project links with user and project joined
        q = (
            db.session.query(
                User.id.label("user_id"),
                User.email,
                User.first_name,
                User.last_name,
                User.auth_id,
                User.is_active,
                User.is_admin,
                Project.id.label("project_id"),
                Project.name.label("project_name"),
                Project.path.label("project_path"),
                Project.owner.label("project_owner_id"),
                Project.is_public,
                Project.is_deleted,
                Project.description.label("project_description"),
                UserProject.can_read,
                UserProject.can_write,
                UserProject.is_owner,
            )
            .join(UserProject, User.id == UserProject.user_id)
            .join(Project, Project.id == UserProject.project_id)
            .order_by(User.email, Project.name)
        )
        rows = q.all()
        auth0_login_stats = _fetch_auth0_login_stats() if include_auth0_login_stats else {}

    columns = [
        "user_id", "email", "first_name", "last_name", "auth_id",
        "is_active", "is_admin",
        "project_id", "project_name", "project_path", "project_owner_id",
        "is_public", "is_deleted", "project_description",
        "can_read", "can_write", "is_owner",
    ]
    if include_auth0_login_stats:
        columns.extend(["last_login", "logins_count"])

    if not rows:
        return pd.DataFrame(columns=columns)

    records = []
    for r in rows:
        record = {
            "user_id": r.user_id,
            "email": r.email,
            "first_name": r.first_name,
            "last_name": r.last_name,
            "auth_id": r.auth_id,
            "is_active": r.is_active,
            "is_admin": r.is_admin,
            "project_id": r.project_id,
            "project_name": r.project_name,
            "project_path": r.project_path,
            "project_owner_id": r.project_owner_id,
            "is_public": r.is_public,
            "is_deleted": r.is_deleted,
            "project_description": r.project_description,
            "can_read": r.can_read,
            "can_write": r.can_write,
            "is_owner": r.is_owner,
        }
        if include_auth0_login_stats:
            stats = auth0_login_stats.get(r.auth_id, {})
            record["last_login"] = stats.get("last_login")
            record["logins_count"] = stats.get("logins_count")
        records.append(record)

    return pd.DataFrame(records, columns=columns)


def main():
    parser = argparse.ArgumentParser(
        description="Collect users and projects from the database and merge into a DataFrame."
    )
    parser.add_argument(
        "--csv",
        metavar="FILE",
        help="Save the merged DataFrame to a CSV file",
    )
    parser.add_argument(
        "--no-print",
        action="store_true",
        help="Do not print the DataFrame to stdout",
    )
    parser.add_argument(
        "--with-auth0-login-stats",
        action="store_true",
        help="Include Auth0 last_login and logins_count fields",
    )
    args = parser.parse_args()

    df = load_users_projects_df(include_auth0_login_stats=args.with_auth0_login_stats)

    if not args.no_print:
        pd.set_option("display.max_columns", None)
        pd.set_option("display.width", None)
        pd.set_option("display.max_colwidth", 50)
        print(df.to_string())

    if args.csv:
        path = Path(args.csv)
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False)
        print(f"\nSaved {len(df)} rows to {path}", file=sys.stderr)

    return df


if __name__ == "__main__":
    main()
