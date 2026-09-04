#!/usr/bin/env bash
#
# deploy_bia.sh - build the sqlite image locally and put it on the bia instances.
#
# Runs from a local machine over one multiplexed SSH connection, so you type
# your password once.
#
# Why this exists: the manual sequence fails silently in two ways.
#
# The tag moves but the containers stay on the old image. `docker compose
# images` then shows <none> against a stale ID and nothing reports an error, so
# a deploy that changed nothing looks like a deploy that worked. This script
# refuses to report success unless every container is running the image ID it
# just loaded.
#
# An arm64 image takes all three instances down with an exec format error. The
# architecture is checked here before anything is transferred, and again on bia
# after the load.
#
# Not to be confused with ./deploy.sh, which does named deployments with
# generated compose overrides and is unrelated to bia.

set -euo pipefail

COMPOSE_FILE=${MDV_COMPOSE_FILE:-/bia/admin_space/sqlite/docker-compose.sqlite-multi.yml}
DEFAULT_SERVICES="mdv_sqlite_1 mdv_sqlite_2 mdv_sqlite_3"
IMAGE=mdvadmin/mdv:sqlite
PLATFORM=linux/amd64
ARCH=amd64
TARBALL=mdv-sqlite-amd64.tar.gz

prog=$(basename "$0")
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

die() { printf '%s: error: %s\n' "$prog" "$*" >&2; exit 1; }
info() { printf '==> %s\n' "$*"; }
warn() { printf 'warning: %s\n' "$*" >&2; }

usage() {
    cat <<EOF
Usage: $prog [-B] [-k] <user@host> [service ...]

Build $IMAGE for $PLATFORM, transfer it to bia, load it, recreate the named
services, and prove the containers actually picked up the new image.

Arguments:
  user@host   ssh target, for example mshaikh@bia.cmd.ox.ac.uk
  service     services to recreate (default: $DEFAULT_SERVICES)

  Only the named services are recreated, so mdv_public is never swept up.

Options:
  -B  skip the build and reuse the $IMAGE already in the local docker
  -k  keep the image tarball on both machines instead of deleting it
  -h  show this help

Environment:
  MDV_COMPOSE_FILE  compose file on the server
                    (default: $COMPOSE_FILE)

Examples:
  $prog mshaikh@bia.cmd.ox.ac.uk
  $prog mshaikh@bia.cmd.ox.ac.uk mdv_sqlite_1
  $prog -B mshaikh@bia.cmd.ox.ac.uk mdv_sqlite_1 mdv_sqlite_2
EOF
}

skip_build=0
keep=0
while getopts ':Bkh' opt; do
    case "$opt" in
        B) skip_build=1 ;;
        k) keep=1 ;;
        h) usage; exit 0 ;;
        :) die "option -$OPTARG needs an argument" ;;
        \?) die "unknown option -$OPTARG (try -h)" ;;
    esac
done
shift $((OPTIND - 1))

[ $# -ge 1 ] || { usage >&2; exit 2; }

target=$1
shift
services=${*:-$DEFAULT_SERVICES}

case "$target" in
    ?*@?*) ;;
    *) die "target must look like user@host, got '$target'" ;;
esac

for svc in $services; do
    case "$svc" in
        mdv_public) die "refusing to recreate mdv_public" ;;
        *[!A-Za-z0-9._-]*) die "invalid service name: '$svc'" ;;
    esac
done

command -v docker >/dev/null || die "docker is not installed"
docker info >/dev/null 2>&1 || die "the docker daemon is not running"
command -v rsync >/dev/null || die "rsync is not installed"

tarball="${TMPDIR:-/tmp}/$TARBALL"
sock=$(mktemp -u "${TMPDIR:-/tmp}/mdv-ssh-XXXXXX")

cleanup() {
    if [ "$keep" -eq 0 ] && [ -f "$tarball" ]; then
        rm -f "$tarball"
    fi
    [ -S "$sock" ] && ssh -S "$sock" -O exit "$target" >/dev/null 2>&1
    return 0
}
trap cleanup EXIT

# --- build -----------------------------------------------------------------

if [ "$skip_build" -eq 0 ]; then
    info "building $IMAGE for $PLATFORM (this takes a while)"
    docker build --platform "$PLATFORM" -t "$IMAGE" "$here" \
        || die "the build failed"
else
    info "skipping the build, using the local $IMAGE"
    docker image inspect "$IMAGE" >/dev/null 2>&1 \
        || die "$IMAGE is not in the local docker, so -B has nothing to reuse"
fi

# An arm64 image starts, then dies with an exec format error, and it takes
# every instance down with it. Catch it before spending time on a 2GB transfer.
built_arch=$(docker image inspect --format '{{.Architecture}}' "$IMAGE")
[ "$built_arch" = "$ARCH" ] || die "$IMAGE is $built_arch, not $ARCH.
Rebuild it with:  docker build --platform $PLATFORM -t $IMAGE ."

built_id=$(docker image inspect --format '{{.Id}}' "$IMAGE")
info "built $ARCH image ${built_id#sha256:}"

# --- save ------------------------------------------------------------------

info "saving and compressing to $tarball"
docker save "$IMAGE" | gzip > "$tarball" || die "docker save failed"

local_sum=$(shasum -a 256 "$tarball" | cut -d' ' -f1)
info "$(du -h "$tarball" | cut -f1), sha256 ${local_sum:0:16}..."

# --- connect ---------------------------------------------------------------

info "connecting to $target"
ssh -f -N -M -S "$sock" -o ControlPersist=60 "$target" \
    || die "could not open an ssh connection to $target"

# ssh joins the command arguments it is given into one string and hands that
# string to the remote login shell, so anything left unquoted here is shell
# syntax on bia rather than data. MDV_COMPOSE_FILE comes from the environment
# and reaches several remote blocks, so it is quoted like everything else.
shquote() {
    local arg quoted=""
    for arg in "$@"; do
        quoted="$quoted '${arg//\'/\'\\\'\'}'"
    done
    printf '%s' "${quoted# }"
}

# Takes a command and its arguments, never a shell snippet. Anything needing a
# pipeline or a redirection goes through an explicit sh -c with its own quoting.
remote() { ssh -S "$sock" "$target" "$(shquote "$@")"; }

remote sh -c 'command -v docker >/dev/null' \
    || die "docker is not available on $target"
remote test -f "$COMPOSE_FILE" \
    || die "no compose file at $COMPOSE_FILE on $target
Point somewhere else with:  MDV_COMPOSE_FILE=... $prog $target"

# --- transfer --------------------------------------------------------------

# --partial --inplace so an interrupted 2GB transfer resumes rather than
# starting over. The tarball lands in the home directory; docker load does not
# care where the file sits, and the compose directory may not be writable.
info "transferring to $target:~/$TARBALL"
rsync -avP --partial --inplace -e "ssh -S $sock" "$tarball" "$target:$TARBALL" \
    || die "the transfer failed"

remote_sum=$(remote sh -c 'sha256sum "$1" | cut -d" " -f1' sh "$TARBALL")
if [ "$local_sum" != "$remote_sum" ]; then
    die "checksum mismatch after transfer.
  local   $local_sum
  remote  $remote_sum
Delete ~/$TARBALL on $target and run this again."
fi
info "checksum matches"

# --- record what the containers run now ------------------------------------

info "recording the current image of each service"
before=$(remote bash -s -- "$COMPOSE_FILE" $services <<'REMOTE'
set -euo pipefail
compose=$1; shift
cd "$(dirname "$compose")"
for svc in "$@"; do
    # || true because compose exits non-zero on an unknown service, and under
    # set -e that would kill the block before the fallback below runs.
    cid=$(docker compose -f "$compose" ps -q "$svc" 2>/dev/null | head -1 || true)
    # Fall back to the container name, which is how the bia compose file names
    # these, in case the compose project name is not derived from the directory.
    [ -n "$cid" ] || cid=$(docker ps --filter "name=^${svc}$" --format '{{.ID}}' | head -1)
    if [ -n "$cid" ]; then
        echo "$svc=$(docker inspect -f '{{.Image}}' "$cid")"
    else
        echo "$svc="
    fi
done
REMOTE
) || die "could not inspect the running services on $target"

printf '%s\n' "$before" | sed 's/^/    /'

# --- load ------------------------------------------------------------------

info "loading the image on $target"
loaded_id=$(remote bash -s -- "$TARBALL" "$IMAGE" "$ARCH" <<'REMOTE'
set -euo pipefail
tarball=$1; image=$2; want_arch=$3

gunzip -c "$tarball" | docker load >&2

arch=$(docker image inspect --format '{{.Architecture}}' "$image")
if [ "$arch" != "$want_arch" ]; then
    echo "loaded image is $arch, not $want_arch" >&2
    exit 1
fi
docker image inspect --format '{{.Id}}' "$image"
REMOTE
) || die "docker load failed on $target"

info "loaded ${loaded_id#sha256:}"

if [ "$loaded_id" != "$built_id" ]; then
    die "the image now tagged $IMAGE on $target is not the one built here.
  built   $built_id
  loaded  $loaded_id
docker load preserves the image ID, so these differ only if something retagged
$IMAGE on $target between the load and this check. The final check compares the
containers against the loaded ID, so carrying on would deploy that other image
and then call it a success. Nothing has been recreated. Find out what moved the
tag, then run this again."
fi

# --- recreate --------------------------------------------------------------

# The compose file sets pull_policy: never, so the load above is what makes
# this possible. Only the named services are passed, so nothing else restarts.
info "recreating: $services"
remote bash -s -- "$COMPOSE_FILE" $services <<'REMOTE'
set -euo pipefail
compose=$1; shift
cd "$(dirname "$compose")"
docker compose -f "$compose" up -d --force-recreate "$@"
REMOTE

# --- prove the containers actually moved -----------------------------------

info "checking that the containers picked up the new image"
after=$(remote bash -s -- "$COMPOSE_FILE" $services <<'REMOTE'
set -euo pipefail
compose=$1; shift
cd "$(dirname "$compose")"
for svc in "$@"; do
    # || true because compose exits non-zero on an unknown service, and under
    # set -e that would kill the block before the fallback below runs.
    cid=$(docker compose -f "$compose" ps -q "$svc" 2>/dev/null | head -1 || true)
    # Fall back to the container name, which is how the bia compose file names
    # these, in case the compose project name is not derived from the directory.
    [ -n "$cid" ] || cid=$(docker ps --filter "name=^${svc}$" --format '{{.ID}}' | head -1)
    if [ -n "$cid" ]; then
        echo "$svc=$(docker inspect -f '{{.Image}}' "$cid")"
    else
        echo "$svc="
    fi
done
REMOTE
) || die "could not inspect the running services on $target"

stale=""
for svc in $services; do
    running=$(printf '%s\n' "$after" | sed -n "s/^${svc}=//p")
    was=$(printf '%s\n' "$before" | sed -n "s/^${svc}=//p")
    if [ -z "$running" ]; then
        stale="$stale $svc(not running)"
    elif [ "$running" != "$loaded_id" ]; then
        stale="$stale $svc"
    fi
    printf '    %-16s %s -> %s\n' "$svc" "${was:-none}" "${running:-none}" \
        | sed 's/sha256://g'
done

if [ -n "$stale" ]; then
    die "these services are not running the image just loaded:$stale
The tag moved but the containers did not. Check for a failed start with:
  ssh $target docker compose -f $COMPOSE_FILE ps"
fi
info "every service is on ${loaded_id#sha256:}"

# --- check the startup sync ------------------------------------------------

# RUN_AUTH0_SYNC_ON_START going missing is silent: the container comes up fine
# and the caches quietly go stale.
info "checking the Auth0 startup sync"
sync_report=$(remote bash -s -- "$COMPOSE_FILE" $services <<'REMOTE'
set -euo pipefail
compose=$1; shift
cd "$(dirname "$compose")"
for svc in "$@"; do
    # || true because compose exits non-zero on an unknown service, and under
    # set -e that would kill the block before the fallback below runs.
    cid=$(docker compose -f "$compose" ps -q "$svc" 2>/dev/null | head -1 || true)
    # Fall back to the container name, which is how the bia compose file names
    # these, in case the compose project name is not derived from the directory.
    [ -n "$cid" ] || cid=$(docker ps --filter "name=^${svc}$" --format '{{.ID}}' | head -1)
    [ -n "$cid" ] || { echo "$svc: not running"; continue; }
    flag=$(docker exec "$cid" printenv RUN_AUTH0_SYNC_ON_START 2>/dev/null || true)
    line=$(docker logs "$cid" 2>&1 \
        | grep -iE "sync on startup|syncing users|sync complete|sync failed" \
        | tail -1 || true)
    echo "$svc: RUN_AUTH0_SYNC_ON_START=${flag:-unset} ${line:+| $line}"
done
REMOTE
) || warn "could not read the sync state"

printf '%s\n' "$sync_report" | sed 's/^/    /'

case "$sync_report" in
    *"RUN_AUTH0_SYNC_ON_START=unset"*)
        warn "the startup sync flag is not set on at least one service.
It lives in the .env.sqlite_N files next to $COMPOSE_FILE. Without it the user
caches go stale and new Auth0 users never appear." ;;
esac

# --- tidy up ---------------------------------------------------------------

if [ "$keep" -eq 0 ]; then
    remote rm -f "$TARBALL" || warn "could not remove ~/$TARBALL on $target"
else
    info "kept $tarball and $target:~/$TARBALL"
fi

printf '\n%s deployed to %s on %s.\n' "$IMAGE" "$services" "$target"
