#!/usr/bin/env bash
#
# onboard.sh - add a user to a bia sqlite MDV instance.
#
# Runs from a local machine. The containers run on bia, so every container step
# goes over one multiplexed SSH connection and you type your password once.
#
# Why this exists: the tenant and connection are read out of the running
# container rather than from a local config.
# The three instances are independent. Run this once per instance.

set -euo pipefail

COMPOSE_SERVICES="mdv_sqlite_1 mdv_sqlite_2 mdv_sqlite_3"
DB_PATH_DEFAULT=/app/mdv/mdv.sqlite3
APP_DIR=/app/python
PERMS_SCRIPT=mdvtools/scripts/manage_project_permissions.py

prog=$(basename "$0")
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

die() { printf '%s: error: %s\n' "$prog" "$*" >&2; exit 1; }
info() { printf '==> %s\n' "$*"; }
warn() { printf 'warning: %s\n' "$*" >&2; }

usage() {
    cat <<EOF
Usage: $prog [-y] <user@host> <service> <email> [project] [permission]

Add a user to a bia sqlite MDV instance: create the Auth0 account, sync it into
the instance database, optionally grant project access, then read the rows back.

Arguments:
  user@host    ssh target, for example mshaikh@bia.cmd.ox.ac.uk
  service      instance to onboard onto ($COMPOSE_SERVICES)
  email        the person's email address
  project      project to grant access to; omitted means sync only
  permission   view, edit or owner (default: view)

Options:
  -y  do not ask for confirmation before touching Auth0
  -h  show this help

Environment:
  M2M_CLIENT_ID      Auth0 machine-to-machine client id
  M2M_CLIENT_SECRET  its secret

  The M2M application needs read:users and create:users on the Management API.
  Credentials come from the environment so they are never written to the repo.

Examples:
  $prog mshaikh@bia.cmd.ox.ac.uk mdv_sqlite_1 furquan.nawab@well.ox.ac.uk
  $prog mshaikh@bia.cmd.ox.ac.uk mdv_sqlite_1 a.person@well.ox.ac.uk cg edit
EOF
}

assume_yes=0
while getopts ':yh' opt; do
    case "$opt" in
        y) assume_yes=1 ;;
        h) usage; exit 0 ;;
        :) die "option -$OPTARG needs an argument" ;;
        \?) die "unknown option -$OPTARG (try -h)" ;;
    esac
done
shift $((OPTIND - 1))

if [ $# -lt 3 ] || [ $# -gt 5 ]; then
    usage >&2
    exit 2
fi

target=$1
service=$2
email=$3
project=${4:-}
permission=${5:-view}

# --- validate the arguments ------------------------------------------------

case "$target" in
    ?*@?*) ;;
    *) die "target must look like user@host, got '$target'" ;;
esac

case " $COMPOSE_SERVICES " in
    *" $service "*) ;;
    *) die "unknown service '$service'; expected one of: $COMPOSE_SERVICES" ;;
esac

case "$email" in
    ?*@?*.?*) ;;
    *) die "'$email' does not look like an email address" ;;
esac

case "$permission" in
    view|edit|owner) ;;
    *) die "permission must be view, edit or owner, got '$permission'" ;;
esac

if [ -n "$project" ]; then
    case "$project" in
        *[!A-Za-z0-9._-]*) die "invalid project name: '$project'" ;;
    esac
fi

# `$prog host service someone@x.ac.uk owner` reads owner as the project name and
# leaves the permission at its default, which is not what anyone means by it.
if [ $# -lt 5 ] && [ -n "$project" ]; then
    case "$project" in
        view|edit|owner)
            die "'$project' is a permission, but it is in the project position.
The order is:  $prog <user@host> <service> <email> [project] [permission]
You may mean:  $prog $target $service $email <project> $project
If a project really is called '$project', name the permission as well." ;;
    esac
fi

# Credentials come from the environment so they are never written to the repo.
[ -n "${M2M_CLIENT_ID:-}" ] || die "set M2M_CLIENT_ID, the Auth0 machine-to-machine client id"
[ -n "${M2M_CLIENT_SECRET:-}" ] || die "set M2M_CLIENT_SECRET, the machine-to-machine secret"

# --- find the local pieces we need -----------------------------------------

# auth_invitation only needs `requests` and the Auth0 Management API over
# HTTPS, so it runs here rather than on bia.
invitation="$here/python/mdvtools/auth/auth_invitation"
[ -f "$invitation" ] || die "cannot find auth_invitation at $invitation
Run this script from its place in the MDV repo."

if [ -x "$here/python/.venv/bin/python" ]; then
    python="$here/python/.venv/bin/python"
elif command -v python3 >/dev/null 2>&1; then
    python=python3
else
    die "no python3 on PATH and no venv at python/.venv"
fi

"$python" -c 'import requests' 2>/dev/null \
    || die "$python cannot import requests, which auth_invitation needs.
Install it with:  $python -m pip install requests"

# --- open one ssh connection -----------------------------------------------

# A single multiplexed connection means the password is typed once rather than
# once per command. BatchMode is deliberately left off so password auth works.
sock=$(mktemp -u "${TMPDIR:-/tmp}/mdv-ssh-XXXXXX")
cfg=""

cleanup() {
    [ -n "$cfg" ] && rm -f "$cfg"
    [ -S "$sock" ] && ssh -S "$sock" -O exit "$target" >/dev/null 2>&1
    return 0
}
trap cleanup EXIT

info "connecting to $target"
ssh -f -N -M -S "$sock" -o ControlPersist=60 "$target" \
    || die "could not open an ssh connection to $target"

# ssh joins the command arguments it is given into one string and hands that
# string to the remote login shell, so anything left unquoted here is shell
# syntax on bia rather than data. Quoting each argument means an email address
# or a service name reaches the remote block as the value it is.
shquote() {
    local arg quoted=""
    for arg in "$@"; do
        quoted="$quoted '${arg//\'/\'\\\'\'}'"
    done
    printf '%s' "${quoted# }"
}

remote() { ssh -S "$sock" "$target" "$(shquote "$@")"; }

# --- read the tenant and connection out of the container -------------------

info "reading the Auth0 tenant and connection from $service"

container_env=$(remote bash -s -- "$service" "$APP_DIR" <<'REMOTE'
set -euo pipefail
service=$1; app_dir=$2

id=$(docker ps --filter "name=^${service}$" --format '{{.ID}}' | head -1)
if [ -z "$id" ]; then
    echo "no running container named $service" >&2
    docker ps --format '  {{.Names}}' >&2
    exit 1
fi
echo "ID=$id"

# printenv first, then the application's own config.json, which is what
# mdv_server_app.load_config does: os.getenv(NAME) or config.get(NAME) against
# the top level of that file. Only the two Auth0 settings get that fallback.
# The app reads MDV_API_ROOT and SQLITE_DB_PATH from the environment alone, so
# a value found in a file would not be the value the app is running with, and
# reporting it would defeat the point of reading the tenant off the container.
for var in AUTH0_DOMAIN AUTH0_DB_CONNECTION MDV_API_ROOT SQLITE_DB_PATH; do
    value=$(docker exec "$id" printenv "$var" 2>/dev/null || true)
    if [ -z "$value" ]; then
        case "$var" in
            AUTH0_DOMAIN|AUTH0_DB_CONNECTION)
                value=$(docker exec "$id" python3 -c '
import json, sys
try:
    cfg = json.load(open(sys.argv[2]))
except Exception:
    sys.exit(0)
found = cfg.get(sys.argv[1])
if isinstance(found, str):
    print(found)
' "$var" "$app_dir/mdvtools/dbutils/config.json" </dev/null 2>/dev/null || true) ;;
        esac
    fi
    echo "${var}=${value}"
done
REMOTE
) || die "could not read the container environment on $target"

container_id=$(printf '%s\n' "$container_env" | sed -n 's/^ID=//p')
domain=$(printf '%s\n' "$container_env" | sed -n 's/^AUTH0_DOMAIN=//p')
connection=$(printf '%s\n' "$container_env" | sed -n 's/^AUTH0_DB_CONNECTION=//p')
api_root=$(printf '%s\n' "$container_env" | sed -n 's/^MDV_API_ROOT=//p')
db_path=$(printf '%s\n' "$container_env" | sed -n 's/^SQLITE_DB_PATH=//p')
db_path=${db_path:-$DB_PATH_DEFAULT}

# Both come from the environment or the top level of the app's config.json, and
# nowhere else. A value nested inside either would not be one the app reads, so
# an empty result here means the setting is genuinely absent.
config_json=$APP_DIR/mdvtools/dbutils/config.json
[ -n "$domain" ] || die "$service does not expose AUTH0_DOMAIN.
Check it by hand:
  ssh $target docker exec $container_id printenv | grep -i auth0
  ssh $target docker exec $container_id cat $config_json"
[ -n "$connection" ] || die "$service does not expose AUTH0_DB_CONNECTION.
An account created in the wrong connection never reaches this instance, so this
script will not guess. Check it by hand:
  ssh $target docker exec $container_id printenv | grep -i auth0
  ssh $target docker exec $container_id cat $config_json"

# --- confirm before touching Auth0 -----------------------------------------

if [ -n "$project" ]; then
    access_desc="$project ($permission)"
else
    access_desc="sync only, no project access"
fi

cat <<EOF

  instance    $service on $target
  container   $container_id
  tenant      $domain
  connection  $connection
  user        $email
  access      $access_desc

EOF

if [ "$assume_yes" -eq 0 ] && [ -t 0 ]; then
    printf 'Create this user in Auth0? [y/N] '
    read -r reply </dev/tty || reply=""
    case "$reply" in
        [yY]|[yY][eE][sS]) ;;
        *) die "cancelled" ;;
    esac
fi

# --- create the Auth0 account ----------------------------------------------

# mktemp creates the file 0600, so the secret is never world-readable, and the
# EXIT trap removes it whichever way the script ends.
cfg=$(mktemp)
cat > "$cfg" <<EOF
domain=$domain
client_id=$M2M_CLIENT_ID
client_secret=$M2M_CLIENT_SECRET
connection=$connection

$email
EOF

info "creating $email in connection '$connection'"

# auth_invitation catches per-user errors, prints them and still exits 0, so
# its exit code is not evidence that the account exists. Check the output too,
# and treat the database read-back at the end as the real proof.
invitation_out=$("$python" "$invitation" "$cfg" 2>&1) || {
    printf '%s\n' "$invitation_out" >&2
    die "auth_invitation failed"
}
printf '%s\n' "$invitation_out"

case "$invitation_out" in
    *"Failed to process"*)
        die "auth_invitation could not create $email, see the message above" ;;
esac

# auth_invitation prints the new Auth0 id when it creates an account, and prints
# nothing identifying when one already exists in this connection. Where we have
# the id, the read-back below matches on it rather than on the email alone, so a
# stale row for the same address cannot stand in for the identity just created.
auth_id=$(printf '%s\n' "$invitation_out" \
    | sed -n 's/^Created user .* with id //p' | tail -1)

# --- sync the user into the instance, and grant access ---------------------

if [ -n "$project" ]; then
    info "granting $permission on '$project' to $email"
else
    info "syncing users from Auth0 into $service"
fi

set +e
grant_out=$(remote bash -s -- \
    "$container_id" "$APP_DIR" "$PERMS_SCRIPT" "$email" "$project" "$permission" <<'REMOTE' 2>&1
set -euo pipefail
id=$1; app_dir=$2; script=$3; email=$4; project=$5; permission=$6

# The image installs its dependencies into a uv-managed venv under $app_dir and
# starts the app with `uv run`. Plain `python` on PATH is the base image's
# interpreter, which cannot import h5py, so mdvtools fails to import before
# argparse ever sees the subcommand. Pick an interpreter that can import the
# package rather than assuming which one that is.
# No -i, and stdin from /dev/null. These blocks reach bia as the script of a
# `bash -s`, which reads them from stdin, and `docker exec -i` drains whatever
# stdin still holds. An -i here swallows the rest of this block, so the commands
# below never run and the block exits 0 having done nothing.
python_bin=$(docker exec "$id" sh -c '
for py in "$0/.venv/bin/python" "$0/.venv/bin/python3" python3 python; do
    if "$py" -c "import mdvtools" >/dev/null 2>&1; then
        printf "%s" "$py"
        exit 0
    fi
done
exit 1' "$app_dir" </dev/null) || {
    echo "nothing in $id can import mdvtools." >&2
    echo "Looked for $app_dir/.venv/bin/python, then python3 and python on PATH." >&2
    echo "List what is there with:  docker exec $id ls $app_dir/.venv/bin" >&2
    exit 1
}

if [ -n "$project" ]; then
    docker exec "$id" sh -c \
        'cd "$0" && "$1" "$2" assign --email "$3" --project "$4" --permission "$5"' \
        "$app_dir" "$python_bin" "$script" "$email" "$project" "$permission" </dev/null
else
    docker exec "$id" sh -c \
        'cd "$0" && "$1" "$2" sync' \
        "$app_dir" "$python_bin" "$script" </dev/null
fi
REMOTE
)
grant_status=$?
set -e

printf '%s\n' "$grant_out"

if [ "$grant_status" -ne 0 ]; then
    case "$grant_out" in
        *"invalid choice: 'sync'"*)
            die "the image running $service predates the 'sync' subcommand of
$PERMS_SCRIPT, which exists on the sqlite branches but not on main.
Rebuild and redeploy the image with ./deploy_bia.sh $target $service, or pass a
project so the 'assign' path is used instead." ;;
    esac
    die "granting access failed on $service (exit $grant_status)"
fi

# --- read the rows back ----------------------------------------------------

info "reading the database rows back"

rows=$(remote bash -s -- "$container_id" "$email" "$db_path" "$auth_id" <<'REMOTE'
set -euo pipefail
id=$1; email=$2; db_path=$3; auth_id=$4

docker exec "$id" python3 -c '
import sqlite3, sys
email, db_path, want_auth_id = sys.argv[1], sys.argv[2], sys.argv[3]
c = sqlite3.connect(db_path)
users = c.execute(
    "select id,email,auth_id,is_admin from users where email=?", (email,)
).fetchall()
projects = c.execute(
    "select user_id,project_id,can_read,can_write,is_owner from user_projects "
    "where user_id in (select id from users where email=?)", (email,)
).fetchall()
print("USERS=%d" % len(users))
for row in users:
    print("  user:", row)
if want_auth_id:
    print("MATCHED=%d" % sum(1 for row in users if row[2] == want_auth_id))
print("PROJECTS=%d" % len(projects))
for row in projects:
    print("  project:", row)
' "$email" "$db_path" "$auth_id" </dev/null
REMOTE
) || die "could not read the database on $service"

printf '%s\n' "$rows"

user_count=$(printf '%s\n' "$rows" | sed -n 's/^USERS=//p')
project_count=$(printf '%s\n' "$rows" | sed -n 's/^PROJECTS=//p')
matched_count=$(printf '%s\n' "$rows" | sed -n 's/^MATCHED=//p')

if [ "${user_count:-0}" -eq 0 ]; then
    die "$email is not in the users table on $service.
The account exists in Auth0 but the sync did not pick it up, which usually means
it landed in a different connection. Compare the user's Identities in the Auth0
dashboard against '$connection'."
fi

if [ -n "$auth_id" ] && [ "${matched_count:-0}" -eq 0 ]; then
    die "$email is in the users table on $service, but no row carries the Auth0
id just created ($auth_id). The rows above belong to an older identity, and
login matches on auth_id rather than email, so this person still cannot sign in.
Compare their Identities in the Auth0 dashboard against '$connection'."
fi

if [ -n "$project" ] && [ "${project_count:-0}" -eq 0 ]; then
    die "$email has no user_projects rows, so the catalog will be empty."
fi

if [ -z "$project" ] && [ "${project_count:-0}" -eq 0 ]; then
    warn "$email has no project access yet, so their catalog will be empty.
Grant some with:  $prog $target $service $email <project> view"
fi

# --- what the operator still has to do -------------------------------------

printf '\n%s is on %s.\n' "$email" "$service"

if [ -n "$api_root" ]; then
    host=${target#*@}
    printf 'Instance URL: https://%s%s\n' "$host" "$api_root"
fi

cat <<EOF

auth_invitation sets a random password and sends nothing, so tell them to reset
it from the sign-in page, or set one in the Auth0 dashboard and pass it on out
of band. Nothing is emailed automatically.
EOF
