from __future__ import annotations

import importlib.util


def _has_module(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


collect_ignore: list[str] = []
collect_ignore_glob: list[str] = []

has_backend_dependencies = all(
    _has_module(module_name)
    for module_name in ("flask_sqlalchemy", "psycopg2")
)
has_auth_dependencies = all(
    _has_module(module_name)
    for module_name in ("authlib", "auth0", "jose", "flask_jwt_extended", "redis")
)

if not has_backend_dependencies:
    collect_ignore.append("mdvtools/dbutils/test")
    collect_ignore_glob.append("mdvtools/dbutils/test/**/*.py")

if not (has_backend_dependencies and has_auth_dependencies):
    collect_ignore.append("mdvtools/auth/tests")
    collect_ignore_glob.append("mdvtools/auth/tests/**/*.py")
