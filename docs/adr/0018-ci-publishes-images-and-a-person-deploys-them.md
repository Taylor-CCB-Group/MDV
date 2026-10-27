# CI publishes images and a person deploys them by digest

Status: proposed

## Context and decision

MDV is a public repository. After each merge to `main`, the "Build and Deploy" workflow pushes `mdvadmin/mdv:<full commit hash>` to Docker Hub and moves `mdvadmin/mdv:latest` to the same image. Each image carries the git build args and OCI labels, so both the running app and `docker inspect` show the commit. A second workflow moves `mdvadmin/mdv:stable` when someone asks for it. The Postgres servers, oceanus and oranus, deploy these images by digest with a script that a person runs on the server.

The SQLite instances on bia, `mdv_sqlite_1` to `mdv_sqlite_3`, are development servers. They run `mdvadmin/mdv:sqlite`. `deploy_bia.sh` builds that image from a local checkout, copies it to bia as a 2 GB tar and loads it, and the compose file sets `pull_policy: never`. Every deploy overwrites the tag, so nobody can redeploy the previous build, and the image does not record which commit it came from. Now that the SQLite backend is on `main`, an image built from `main` uses SQLite whenever `DB_BACKEND` is not set, so previews on bia can run the same images as the Postgres servers.

Nothing builds an image for a pull request, and no workflow tests an image before it reaches a server.

An earlier plan had CI deploy to bia and oceanus itself, through self-hosted runners on each server, with GitHub Environment approvals as the gate.

We decided instead that CI builds every image and publishes it to Docker Hub, and a person deploys it by digest with a script. There are no self-hosted runners.

## How it works

### Images

The Docker Hub workflow keeps building `main` as it does now. A second workflow builds a pull request opened from this repository only when it carries the `preview` label, and runs again on each new commit while the label stays. It first waits for the pull request's normal Python and TypeScript CI runs on that commit, so pytest, pyright, the TypeScript type check, the Vite build and Vitest must all pass before anything is built. It waits for the Python run only when the pull request changes `python/`, because that is when Python CI runs. It then builds the head commit for amd64 and pushes it as `mdvadmin/mdv:pr-<full commit hash>`. The prefix keeps a bare hash tag meaning a commit of `main`, and a pull request build never moves `:latest`. Pull request builds read the shared build cache but do not write to it, so a branch cannot replace the cache that `main` builds reuse. Pull requests from forks get no image, because GitHub does not pass the Docker Hub secrets to their workflows.

Playwright runs against each `main` image inside the CI job, with `ENABLE_AUTH=0`. Publishing a GitHub Release adds a `v<version>` tag to the existing image for that commit if its Playwright run passed. Nothing is rebuilt for a release. The `:stable` workflow stays as it is.

### Deploys

A person deploys to the Postgres servers with a script that runs on the server. It pulls an image by digest, dumps every database with `pg_dumpall`, and tags the image that each service runs now so it can roll back. It points the image name that the compose file already uses at the pulled image, then recreates only the services it is given. For each one, it checks that the container runs the pulled image and waits for the app to answer HTTP. It then reads the startup log, because `safe_mdv_app.py` falls back to a Flask app with no routes when startup fails, and that app answers HTTP too. The log must also report the Postgres backend. The compose files stay as they are. The script lives in the operations notes and not in this repository, because it names production services and paths.

The SQLite instances on bia keep `deploy_bia.sh`, which a person runs from their own machine, and this decision does not change it. On bia, CI images reach the instances only as previews.

Rolling back means pointing the image name at the rollback tag, or deploying the previous digest. To see what a server is running, read the image itself. Pulling by digest records the digest in the image's `RepoDigests`, and the OCI revision label holds the commit.

### Previews

A preview runs a pull request on one of the SQLite instances on bia, with that instance's Auth0 login and the data already in its volume. `scripts/bia_preview_poller.py` runs from cron under a team member's account on bia every 5 minutes and looks for open pull requests that carry the `preview` label. When one has an image for its head commit and an instance is free, the poller copies that instance's SQLite database with the backup API and pulls the image. It then recreates the service on that image through an override file of its own, so the shared compose file stays untouched. Like the Postgres deploy, it checks the running image, the HTTP answer and the startup log, then comments the instance's link on the pull request.

A new commit redeploys the same instance and edits the same comment. When the pull request closes or loses the label, the poller stops the service, puts the database copy back, and recreates the service on `mdvadmin/mdv:sqlite`. A preview that fails to start is torn down the same way, and the poller waits for a new commit before trying again. When all 3 instances are busy, it comments that the pull request is waiting. The poller uses a fine-grained token that can read this repository's pull requests and workflow runs and write comments on its pull requests.

## Why

GitHub's [security guidance for Actions](https://docs.github.com/en/actions/reference/secure-use-reference) says that self-hosted runners "should almost never be used for public repositories", because anyone can open a pull request. It also says environments with required reviewers do not isolate a self-hosted runner, so the approval gate in the earlier plan would not have reduced that risk. A runner on bia would sit beside the instances people use.

Runners also needed write access to the admin directories on both servers. Oceanus runs several sites from one compose file, so a runner there also needed a decision about where it lives. The server administrators own both. A script that a person runs needs neither, and it keeps working when GitHub Actions is down.

Docker Hub already holds what a digest deploy needs. Every `main` build has a tag named after its commit and labels that record it, and the Postgres servers already pull from there. A second registry would put each `main` commit in two places that have to stay in step. A manual run of the workflow can push the same hash tag again with a new digest, so deploys name the digest, which no later push can move.

Swapping the image back does not undo what the new image wrote at startup, such as the ownership it gives administrators over projects that have no owner. That is why the deploy backs up the database first.

The release gate tests the image inside CI because that is how the Playwright suite works. Its project tests generate data on the test runner's own disk and expect the server to see that data through a bind mount, and the suite has no login handling. So the suite cannot test a deployed instance. Testing in CI also means GitHub never needs a network route to bia. A release retags the image that was tested, so what ships is exactly what passed.

## Considered and rejected

- CI deploying to both servers through runners, as in the earlier plan, for the reasons above.
- CI deploying previews to bia automatically while people deploy production. This still needs a runner on bia.
- Publishing to a public GHCR package alongside Docker Hub. The package would sit under the repository's own permissions, but each `main` commit would exist in two registries, and the Postgres servers already pull from Docker Hub.
- An image for every pull request. Most pull requests never get deployed, so building each one wastes CI time and fills Docker Hub.
- Building from a `workflow_run` trigger once the test workflows finish, so the tests run only once. It runs from `main`'s copy of the workflow with secrets, has to ask the API whether every check passed, and needs separate handling for a label added after the tests finish.
- Calling the Python and TypeScript workflows from the image workflow. Every labelled commit would run its tests twice.
- Merging the Python and TypeScript workflows into one, so the image job can wait on their jobs directly. It needs a third-party action for the Python path filter and changes the check names that branch protection uses.
- One deploy script in this repository for every server. It would have to take production service names and paths as configuration, and the repository is public.
- A compose file that reads the digest from a variable. It is cleaner, but it needs write access to the compose directory.
- Building and testing on each Release. The released image would be a new build that nobody had run before, and a failed test would leave a Release with no image.
- Testing a staging URL. Only the route-mocked catalog tests can run against a remote instance.
- Separate preview containers with dev auth and a synthetic project. They need new ports or nginx paths, and a synthetic project that keeps up with the code, while the SQLite instances already have Auth0, addresses and data.
- Changing `deploy_bia.sh` to pull images by digest. It runs from a person's machine, while previews need a job that lives on bia.
- A keyword in the pull request description as the preview trigger. Any author can add a keyword, but only people with triage or write access can add a label.

## Consequences

- Every app service on the Postgres servers must set `DB_BACKEND=postgres`. Without it, an image from `main` starts on a new, empty SQLite database, which is why the Postgres script stops when a startup log reports another backend.
- Database migrations stay manual. The release notes name any migration, and an operator runs it with the app stopped before deploying that release.
- Retagging moves a tag for every service on the server that uses it. A service that the script was not asked to recreate picks up the new image the next time anyone recreates it.
- The Postgres script lives outside the repository, so each person's copy on a server can drift from the operations notes.
- The image job holds a runner while it waits for the tests, for up to an hour.
- Re-running a failed test job does not build an image. Pushing a commit, or removing and adding the label, starts a new image run.
- The tests run on GitHub's merge commit with `main`, while the image is built from the head commit.
- Pull request images stay on Docker Hub, and nothing removes them.
- Docker Hub limits how often anonymous clients can pull. If the preview poller reaches that limit, it can log in with a read-only Docker Hub token.
- `playwright.yml` has never passed, and it has to be fixed before the release gate can work.
- Nobody can count on a SQLite instance on bia running `:sqlite` while a preview is labelled.
- A preview runs unmerged code against an instance's real data. The database is restored afterwards, but edits to project files are not.
- Running `deploy_bia.sh` over an instance that holds a preview puts it on `:sqlite` without restoring its database. The poller's next run deploys the preview again.
