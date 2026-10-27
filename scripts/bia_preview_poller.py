#!/usr/bin/env python3
"""Run pull requests labelled preview on the bia SQLite instances (ADR 0018).

Cron runs this on bia every 5 minutes under a team member's account:

    */5 * * * * python3 ~/mdv-preview/bia_preview_poller.py >> ~/.config/mdv-preview/poller.log 2>&1

Each run asks GitHub for open pull requests into main that carry the preview label. The services named in the config are the preview slots. A labelled pull
request whose head commit has a Docker Hub image (built by pull-request-image.yml once its tests pass) takes a free slot:
1. The poller pulls mdvadmin/mdv:pr-<head commit>.
2. Before the slot's first preview, it copies the instance's SQLite database
   with the backup API, next to the database on the volume.
3. It tags the image mdvadmin/mdv:preview-<service> and recreates the service
   with an override file that points it at that tag. The shared compose file
   is never edited.
4. It checks that the container runs the image, answers HTTP and logs a clean
   start on SQLite, then comments the link on the pull request.

A new commit redeploys the same slot. When the label comes off or the pull
request closes, the poller stops the service, puts the database copy back and
recreates the service on mdvadmin/mdv:sqlite. Edits to project files made
during a preview are not rolled back.

Running deploy_bia.sh over a slot puts it back on :sqlite without restoring the
database. The next run sees that and deploys the preview again.

Config, in ~/.config/mdv-preview/config.json, readable only by its owner:

    {
      "token": "github_pat_...",
      "compose_file": "/bia/admin_space/sqlite/docker-compose.sqlite-multi.yml",
      "slots": {
        "mdv_sqlite_1": "https://<address of instance 1>/",
        "mdv_sqlite_2": "https://<address of instance 2>/",
        "mdv_sqlite_3": "https://<address of instance 3>/"
      }
    }

The token is a fine-grained token for this repository only, with read and
write access to pull requests and read access to actions.

Run with --dry-run to print what the poller would do without changing
anything on bia or GitHub.
"""

import argparse
import fcntl
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request

REPO = "Taylor-CCB-Group/MDV"
LABEL = "preview"
WORKFLOW = "pull-request-image.yml"
IMAGE = "mdvadmin/mdv"
BASE_TAG = IMAGE + ":sqlite"
PLATFORM = "linux/amd64"
ARCH = "amd64"
STARTUP_TIMEOUT = 300
CONTAINER_PYTHON = "/app/python/.venv/bin/python"
DEFAULT_DB = "/app/mdv/mdv.sqlite3"
STATE_DIR = os.path.expanduser("~/.config/mdv-preview")
# Lets the poller find its own comment again if the state file is lost.
MARKER = "<!-- mdv-preview -->"

# safe_mdv_app.py serves a Flask app with no routes when startup fails, and
# that app answers HTTP too, so the log decides whether startup worked.
STARTUP_ERRORS = re.compile(
    r"Error during app initialization|Error importing mdv_app|using fallback"
)
SQLITE_SELECTED = "Database backend selected: sqlite"

# Any status code counts, because the root path may redirect to the login page.
HTTP_PROBE = """
import http.client
c = http.client.HTTPConnection("127.0.0.1", 5055, timeout=5)
c.request("GET", "/")
c.getresponse()
"""

# Takes a consistent copy while the app is running, then checks it.
BACKUP = """
import sqlite3, sys
source = sqlite3.connect(sys.argv[1])
copy = sqlite3.connect(sys.argv[2])
source.backup(copy)
source.close()
if copy.execute("pragma quick_check").fetchone()[0] != "ok":
    sys.exit(sys.argv[2] + " failed quick_check")
copy.close()
"""

# Runs in a one-off container while the service is stopped. The WAL files
# belong to the preview's database, so they go too. The snapshot stays until
# the service is back up, so a teardown that fails halfway can run again.
RESTORE = 'cp -p "$1" "$2" && rm -f "$2-wal" "$2-shm"'


class DeployError(Exception):
    pass


class PullError(DeployError):
    """The image could not be fetched, so nothing on the slot changed."""


def log(message):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), message, flush=True)


# --- GitHub -----------------------------------------------------------------


class GitHub:
    def __init__(self, token, dry_run):
        self.token = token
        self.dry_run = dry_run

    def request(self, method, path, body=None):
        url = path if path.startswith("https://") else "https://api.github.com" + path
        data = None if body is None else json.dumps(body).encode()
        req = urllib.request.Request(url, data=data, method=method)
        req.add_header("Authorization", "Bearer " + self.token)
        req.add_header("Accept", "application/vnd.github+json")
        req.add_header("X-GitHub-Api-Version", "2022-11-28")
        req.add_header("User-Agent", "mdv-bia-preview-poller")
        with urllib.request.urlopen(req, timeout=30) as response:
            link = response.headers.get("Link", "")
            return json.loads(response.read().decode() or "null"), link

    def get_all(self, path, key=None):
        """Follows pagination. key names the list in endpoints that wrap it in an object."""
        items = []
        while path:
            page, link = self.request("GET", path)
            items.extend(page if key is None else page[key])
            match = re.search(r'<([^>]+)>; rel="next"', link)
            path = match.group(1) if match else None
        return items

    def labelled_pull_requests(self):
        pulls = self.get_all(
            "/repos/{}/pulls?state=open&base=main&per_page=100".format(REPO)
        )
        return {
            pull["number"]: pull["head"]["sha"]
            for pull in pulls
            if LABEL in [label["name"] for label in pull["labels"]]
            # Forks get no image, so there is nothing to deploy.
            and (pull["head"]["repo"] or {}).get("full_name") == REPO
        }

    def image_status(self, sha):
        """Says whether pull-request-image.yml has built an image for sha."""
        runs, _ = self.request(
            "GET",
            "/repos/{}/actions/workflows/{}/runs?head_sha={}&per_page=30".format(
                REPO, WORKFLOW, sha
            ),
        )
        runs = runs["workflow_runs"]
        if any(run["status"] != "completed" for run in runs):
            return "pending"
        # Adding another label starts a run whose jobs are all skipped, and
        # that run can still succeed, so only a successful image job counts.
        for run in runs:
            jobs = self.get_all(
                "/repos/{}/actions/runs/{}/jobs?per_page=100".format(REPO, run["id"]), "jobs"
            )
            if any(job["name"] == "image" and job["conclusion"] == "success" for job in jobs):
                return "ready"
        # Failed tests skip the image job, so the run says whether they failed.
        if any(run["conclusion"] == "failure" for run in runs):
            return "failed"
        return "missing"

    def comment(self, state, number, body):
        """Keeps one comment per pull request, edited in place."""
        body = body + "\n\n" + MARKER
        known = state["comments"].get(str(number))
        if known is None and not self.dry_run:
            # The state file was lost or reset, so look for the marker before
            # posting a second comment.
            comments = self.get_all("/repos/{}/issues/{}/comments?per_page=100".format(REPO, number))
            found = next((c for c in comments if MARKER in (c["body"] or "")), None)
            if found:
                known = {"id": found["id"], "body": found["body"]}
                state["comments"][str(number)] = known
        if known and known["body"] == body:
            return
        log("comment on #{}: {}".format(number, body.splitlines()[0]))
        if self.dry_run:
            return
        if known:
            try:
                self.request(
                    "PATCH",
                    "/repos/{}/issues/comments/{}".format(REPO, known["id"]),
                    {"body": body},
                )
                known["body"] = body
                return
            except urllib.error.HTTPError as error:
                # The comment was deleted, or someone else wrote the marker
                # into theirs, so post a new one.
                if error.code not in (403, 404):
                    raise
        created, _ = self.request(
            "POST", "/repos/{}/issues/{}/comments".format(REPO, number), {"body": body}
        )
        state["comments"][str(number)] = {"id": created["id"], "body": body}


# --- Docker -----------------------------------------------------------------


class Docker:
    def __init__(self, compose_file, override_file, dry_run):
        self.compose_file = compose_file
        self.override_file = override_file
        self.dry_run = dry_run

    def run(self, args, check=True, change=False):
        if change and self.dry_run:
            log("would run: " + " ".join(args))
            return ""
        result = subprocess.run(
            args,
            cwd=os.path.dirname(self.compose_file),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
        )
        if check and result.returncode != 0:
            raise DeployError("{} failed: {}".format(" ".join(args), result.stdout.strip()))
        return result.stdout.strip()

    def compose(self, *args, override=False, change=True):
        files = ["-f", self.compose_file]
        if override:
            files += ["-f", self.override_file]
        return self.run(["docker", "compose"] + files + list(args), change=change)

    def container_of(self, service):
        cid = self.run(
            ["docker", "compose", "-f", self.compose_file, "ps", "-q", service],
            check=False,
        )
        if not cid:
            # The bia compose file names containers after their services.
            cid = self.run(
                ["docker", "ps", "--filter", "name=^{}$".format(service), "--format", "{{.ID}}"]
            )
        return cid.splitlines()[0] if cid else ""

    def running_image(self, service):
        cid = self.container_of(service)
        if not cid:
            return ""
        return self.run(["docker", "inspect", "-f", "{{.Image}}", cid])

    def image_id(self, ref):
        return self.run(["docker", "image", "inspect", "-f", "{{.Id}}", ref], check=False)

    def pull(self, sha):
        ref = "{}:pr-{}".format(IMAGE, sha)
        try:
            self.run(["docker", "pull", "--platform", PLATFORM, ref])
        except DeployError as error:
            raise PullError(str(error))
        arch = self.run(["docker", "image", "inspect", "-f", "{{.Architecture}}", ref])
        if arch != ARCH:
            raise DeployError("{} is {}, not {}".format(ref, arch, ARCH))
        digests = self.run(["docker", "image", "inspect", "-f", "{{json .RepoDigests}}", ref])
        digest = next(
            (d.split("@", 1)[1] for d in json.loads(digests) if d.startswith(IMAGE + "@")),
            "unknown",
        )
        return self.image_id(ref), digest

    def database_path(self, cid):
        path = self.run(["docker", "exec", cid, "printenv", "SQLITE_DB_PATH"], check=False)
        return path or DEFAULT_DB

    def write_override(self, slots):
        # JSON is valid YAML, so compose reads this as an override file.
        services = {
            service: {"image": preview_tag(service)}
            for service, entry in slots.items()
            if entry.get("sha")
        }
        if self.dry_run:
            return
        with open(self.override_file, "w") as handle:
            json.dump({"services": services}, handle, indent=2)

    def check_started(self, service, expected_id):
        """Raises DeployError unless the service runs expected_id and started cleanly."""
        cid = self.container_of(service)
        running = self.run(["docker", "inspect", "-f", "{{.Image}}", cid]) if cid else ""
        if running != expected_id:
            raise DeployError("{} is running {}, not {}".format(service, running or "nothing", expected_id))
        # Gunicorn binds the port before its worker builds the app, so the
        # first answer comes once startup has finished.
        deadline = time.time() + STARTUP_TIMEOUT
        probe = ["docker", "exec", cid, CONTAINER_PYTHON, "-c", HTTP_PROBE]
        while subprocess.run(probe, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode != 0:
            if time.time() > deadline:
                raise DeployError("{} did not answer HTTP within {}s".format(service, STARTUP_TIMEOUT))
            time.sleep(5)
        logs = self.run(["docker", "logs", cid])
        if STARTUP_ERRORS.search(logs):
            raise DeployError("{} logged a startup error".format(service))
        if SQLITE_SELECTED not in logs:
            raise DeployError('{} did not report "{}"'.format(service, SQLITE_SELECTED))


def preview_tag(service):
    return "{}:preview-{}".format(IMAGE, service)


# --- deploy and tear down ---------------------------------------------------


def deploy(docker, state, service, number, sha):
    slots = state["slots"]
    if docker.dry_run:
        log("would deploy #{} ({}) to {}".format(number, sha[:7], service))
        return "unknown"
    image_id, digest = docker.pull(sha)

    if service not in slots:
        # First preview on this slot: copy the database while the service
        # still runs :sqlite, and record the slot before changing anything, so
        # a failed run still knows to restore it.
        cid = docker.container_of(service)
        if not cid:
            raise DeployError("{} is not running, so its database cannot be copied".format(service))
        db = docker.database_path(cid)
        snapshot = db + ".preview-snapshot"
        docker.run(["docker", "exec", cid, CONTAINER_PYTHON, "-c", BACKUP, db, snapshot], change=True)
        slots[service] = {"pr": number, "db": db, "snapshot": snapshot}
        save_state(state)
        log("copied the database of {} to {}".format(service, snapshot))

    slots[service].update({"pr": number, "sha": sha, "digest": digest, "image_id": image_id})
    docker.run(["docker", "tag", image_id, preview_tag(service)], change=True)
    docker.write_override(slots)
    save_state(state)
    docker.compose("up", "-d", "--no-deps", "--force-recreate", service, override=True)
    docker.check_started(service, image_id)
    log("#{} ({}) runs on {}".format(number, sha[:7], service))
    return digest


def tear_down(docker, state, service):
    entry = state["slots"][service]
    if docker.dry_run:
        log("would put {} back on {} and restore its database".format(service, BASE_TAG))
        return
    docker.compose("stop", service)
    docker.compose(
        "run", "--rm", "--no-deps", "-T", "--entrypoint", "sh", service,
        "-c", RESTORE, "sh", entry["snapshot"], entry["db"],
    )
    entry["sha"] = None
    docker.write_override(state["slots"])
    docker.compose("up", "-d", "--no-deps", "--force-recreate", service)
    docker.check_started(service, docker.image_id(BASE_TAG))
    docker.run(["docker", "exec", docker.container_of(service), "rm", "-f", entry["snapshot"]], check=False, change=True)
    del state["slots"][service]
    save_state(state)
    docker.run(["docker", "rmi", preview_tag(service)], check=False, change=True)
    log("{} is back on {} with its database restored".format(service, BASE_TAG))


# --- one run ----------------------------------------------------------------


def poll(github, docker, config, state):
    slots = config["slots"]
    labelled = github.labelled_pull_requests()
    failures = 0

    # Free the slots of pull requests that closed or lost the label.
    for service in list(state["slots"]):
        number = state["slots"][service]["pr"]
        if number in labelled:
            continue
        try:
            tear_down(docker, state, service)
            github.comment(
                state, number,
                "Preview removed. `{}` is back on `{}` with its database restored.".format(service, BASE_TAG),
            )
        except (DeployError, urllib.error.URLError) as error:
            failures += 1
            log("could not tear down {}: {}".format(service, error))

    # Oldest pull requests get slots first.
    for number in sorted(labelled):
        sha = labelled[number]
        service = next((s for s, e in state["slots"].items() if e["pr"] == number), None)
        entry = state["slots"].get(service, {})
        current = service and entry.get("sha") == sha
        if current and docker.running_image(service) == entry.get("image_id"):
            continue
        # A commit that failed to start is not retried, so the instance is not
        # switched back and forth every run. A new commit tries again.
        if state["failed"].get(str(number)) == sha:
            continue
        try:
            status = "ready" if current else github.image_status(sha)
            running = (
                " The preview of {} keeps running until then.".format(entry["sha"][:7])
                if entry.get("sha") else ""
            )
            if status == "pending":
                github.comment(state, number, "Waiting for the tests and image build for {}.{}".format(sha[:7], running))
                continue
            if status == "failed":
                github.comment(state, number, "No preview for {}, because its tests or image build failed.{}".format(sha[:7], running))
                continue
            if status == "missing":
                github.comment(
                    state, number,
                    "No image was built for {}. Push a commit, or remove and add the `{}` label again.{}".format(sha[:7], LABEL, running),
                )
                continue
            if service is None:
                service = next((s for s in slots if s not in state["slots"]), None)
                if service is None:
                    github.comment(
                        state, number,
                        "Waiting for a free instance. All {} bia SQLite instances are running other previews.".format(len(slots)),
                    )
                    continue
            try:
                digest = deploy(docker, state, service, number, sha)
            except PullError as error:
                # A rate limit or network error says nothing about the commit,
                # and the slot is untouched, so the next run tries again.
                failures += 1
                log("could not pull the image for #{}: {}".format(number, error))
                continue
            except DeployError as error:
                failures += 1
                log("could not deploy #{} to {}: {}".format(number, service, error))
                state["failed"][str(number)] = sha
                outcome = "`{}` was not changed.".format(service)
                if service in state["slots"]:
                    tear_down(docker, state, service)
                    outcome = "`{}` went back to `{}` with its database restored.".format(service, BASE_TAG)
                github.comment(
                    state, number,
                    "The preview of {} did not start: {}. {}".format(sha[:7], error, outcome),
                )
                continue
            github.comment(
                state, number,
                "Preview of {} is running at {}\n\n"
                "Instance `{}`, image `{}@{}`. It updates on each new commit while the `{}` label stays. "
                "When the label comes off or the pull request closes, the instance goes back to `{}` "
                "and its database is restored. Edits to project files are not rolled back.".format(
                    sha[:7], slots[service], service, IMAGE, digest, LABEL, BASE_TAG
                ),
            )
            state["failed"].pop(str(number), None)
        except (DeployError, urllib.error.URLError) as error:
            failures += 1
            log("could not handle #{}: {}".format(number, error))
    return failures


# --- state and entry point --------------------------------------------------


def state_path():
    return os.path.join(STATE_DIR, "state.json")


def load_state():
    state = {"slots": {}, "comments": {}, "failed": {}}
    try:
        with open(state_path()) as handle:
            state.update(json.load(handle))
    except FileNotFoundError:
        pass
    return state


def save_state(state):
    temp = state_path() + ".tmp"
    with open(temp, "w") as handle:
        json.dump(state, handle, indent=2)
    os.replace(temp, state_path())


def load_config():
    path = os.path.join(STATE_DIR, "config.json")
    if not os.path.exists(path):
        sys.exit("no config at {}. See the top of this script for its format.".format(path))
    if os.stat(path).st_mode & 0o077:
        sys.exit("{} holds a token, so make it readable only by you: chmod 600 {}".format(path, path))
    with open(path) as handle:
        config = json.load(handle)
    for service in config["slots"]:
        # Shares the compose file with the slots but is not a dev server.
        if service == "mdv_public":
            sys.exit("refusing to use mdv_public as a preview slot")
        if not re.fullmatch(r"[A-Za-z0-9._-]+", service):
            sys.exit("invalid service name in config: {!r}".format(service))
    return config


def main():
    parser = argparse.ArgumentParser(description="Run pull requests labelled preview on the bia SQLite instances.")
    parser.add_argument("--dry-run", action="store_true", help="print what would change and change nothing")
    args = parser.parse_args()

    os.makedirs(STATE_DIR, mode=0o700, exist_ok=True)
    # A run waits up to STARTUP_TIMEOUT per slot, so it can outlast the cron
    # interval. The next run then exits instead of competing with it.
    lock = open(os.path.join(STATE_DIR, "lock"), "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        return 0

    config = load_config()
    state = load_state()
    github = GitHub(config["token"], args.dry_run)
    docker = Docker(config["compose_file"], os.path.join(STATE_DIR, "override.yml"), args.dry_run)
    # Saved even when poll raises, so comment IDs and failed commits survive.
    try:
        failures = poll(github, docker, config, state)
    finally:
        if not args.dry_run:
            save_state(state)
    # A non-zero exit makes cron mail the log to the account's owner.
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
