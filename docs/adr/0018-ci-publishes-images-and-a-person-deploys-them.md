# CI publishes images and a person deploys them by digest

Status: proposed

## Context and decision

MDV is a public repository. After each merge to `main`, the "Build and Deploy" workflow pushes `mdvadmin/mdv:<full commit hash>` to Docker Hub and moves `mdvadmin/mdv:latest` to the same image. Each image carries the git build args and OCI labels, so both the running app and `docker inspect` show the commit. A second workflow moves `mdvadmin/mdv:stable` when someone asks for it. The Postgres servers, oceanus and oranus, deploy these images by digest with a script that a person runs on the server.

The SQLite test instances on bia, known as the C3Collide instances, run `mdvadmin/mdv:sqlite`. `deploy_bia.sh` builds that image from a local checkout, copies it to bia as a 2 GB tar and loads it, and the compose file sets `pull_policy: never`. Every deploy overwrites the tag, so nobody can redeploy the previous build, and the image does not record which commit it came from. Now that the SQLite backend is on `main`, an image built from `main` uses SQLite whenever `DB_BACKEND` is not set, so bia can run the same images as the Postgres servers.

Nothing builds an image for a pull request, and no workflow tests an image before it reaches a server.

An earlier plan had CI deploy to bia and oceanus itself, through self-hosted runners on each server, with GitHub Environment approvals as the gate.

We decided instead that CI builds every image and publishes it to Docker Hub, and a person deploys it by digest with a script. There are no self-hosted runners.

## How it works

### Images

The Docker Hub workflow keeps building `main` as it does now. It also builds the head of every pull request opened from this repository, for amd64, and pushes it as `mdvadmin/mdv:pr-<full commit hash>`. The prefix keeps a bare hash tag meaning a commit of `main`, and a pull request build never moves `:latest`. Pull request builds read the shared build cache but do not write to it, so a branch cannot replace the cache that `main` builds reuse. Pull requests from forks get no image, because GitHub does not pass the Docker Hub secrets to their workflows.

Playwright runs against each `main` image inside the CI job, with `ENABLE_AUTH=0`. Publishing a GitHub Release adds a `v<version>` tag to the existing image for that commit if its Playwright run passed. Nothing is rebuilt for a release. The `:stable` workflow stays as it is.

### Deploys

Each kind of server has its own deploy script, and both scripts follow the same steps. The script pulls an image by digest, backs up the database, and tags the image that each service runs now so it can roll back. It points the image name that the compose file already uses at the pulled image, then recreates only the services it is given. For each one, it checks that the container runs the pulled image and waits for the app to answer HTTP. It then reads the startup log, because `safe_mdv_app.py` falls back to a Flask app with no routes when startup fails, and that app answers HTTP too. The compose files stay as they are, `pull_policy: never` included.

On the Postgres servers, the script dumps every database with `pg_dumpall` and checks that each startup log reports the Postgres backend. It lives in the operations notes and not in this repository, because it names production services and paths. On bia, `deploy_bia.sh` pulls the digest instead of building a tar. Before it recreates anything, it copies the SQLite database with SQLite's backup API, run inside the container.

Rolling back means pointing the image name at the rollback tag, or deploying the previous digest. To see what a server is running, read the image itself. Pulling by digest records the digest in the image's `RepoDigests`, and the OCI revision label holds the commit.

### Previews

A preview is a throwaway instance of a pull request on bia. A cron job under a team member's account on bia polls GitHub every 5 minutes for open pull requests that carry the `preview` label. For each one, it runs that pull request's image as a single SQLite container with dev auth and a fixed synthetic project, in one of a few reserved port slots, and comments the link on the pull request. When the pull request gets a new commit, the job redeploys it and edits the same comment. When the pull request closes or loses the label, the job removes the preview. When every slot is busy, it comments that the pull request is waiting. The job uses a fine-grained token that can only write pull request comments on MDV.

## Why

GitHub's [security guidance for Actions](https://docs.github.com/en/actions/reference/secure-use-reference) says that self-hosted runners "should almost never be used for public repositories", because anyone can open a pull request. It also says environments with required reviewers do not isolate a self-hosted runner, so the approval gate in the earlier plan would not have reduced that risk. A runner on bia would sit beside the instances people use.

Runners also needed write access to the admin directories on both servers. Oceanus runs several sites from one compose file, so a runner there also needed a decision about where it lives. The server administrators own both. A script that a person runs needs neither, and it keeps working when GitHub Actions is down.

Docker Hub already holds what a digest deploy needs. Every `main` build has a tag named after its commit and labels that record it, and the Postgres servers already pull from there. A second registry would put each `main` commit in two places that have to stay in step. A manual run of the workflow can push the same hash tag again with a new digest, so deploys name the digest, which no later push can move. Pulling by digest also replaces the tar transfer to bia.

The two kinds of server back up different databases, and each startup log must report a different backend. The Postgres script also names production services that this public repository should not list. So each kind of server keeps its own script, and this decision fixes the steps the scripts share.

Swapping the image back does not undo what the new image wrote at startup, such as the ownership it gives administrators over projects that have no owner. That is why every deploy backs up the database first, on bia as well as on the Postgres servers.

The release gate tests the image inside CI because that is how the Playwright suite works. Its project tests generate data on the test runner's own disk and expect the server to see that data through a bind mount, and the suite has no login handling. So the suite cannot test a deployed instance. Testing in CI also means GitHub never needs a network route to bia. A release retags the image that was tested, so what ships is exactly what passed.

## Considered and rejected

- CI deploying to both servers through runners, as in the earlier plan, for the reasons above.
- CI deploying previews to bia automatically while people deploy production. This still needs a runner on bia.
- Publishing to a public GHCR package alongside Docker Hub. The package would sit under the repository's own permissions, but each `main` commit would exist in two registries, and the Postgres servers already pull from Docker Hub.
- One deploy script in this repository for every server. It would have to take production service names and paths as configuration, and the repository is public.
- Skipping the backup on bia because its instances are for testing. The startup writes described above happen there too.
- A compose file that reads the digest from a variable. It is cleaner, but it needs write access to the compose directory.
- Building and testing on each Release. The released image would be a new build that nobody had run before, and a failed test would leave a Release with no image.
- Testing a staging URL. Only the route-mocked catalog tests can run against a remote instance.
- An `auth=auth0` option for previews. Previews use dev auth and synthetic data only, so they need no Auth0 callback slots.
- A keyword in the pull request description as the preview trigger. Any author can add a keyword, but only people with triage or write access can add a label.

## Consequences

- Every app service on the Postgres servers must set `DB_BACKEND=postgres`. Without it, an image from `main` starts on a new, empty SQLite database, which is why the Postgres script stops when a startup log reports another backend.
- Database migrations stay manual. The release notes name any migration, and an operator runs it with the app stopped before deploying that release.
- Retagging moves a tag for every service on the server that uses it. A service that the script was not asked to recreate picks up the new image the next time anyone recreates it.
- The Postgres script lives outside the repository, so each person's copy on a server can drift from the operations notes, and a change to the shared steps has to be made in both scripts.
- Pull request images stay on Docker Hub, and nothing removes them.
- Docker Hub limits how often anonymous clients can pull. If the preview poller reaches that limit, it can log in with a read-only Docker Hub token.
- `playwright.yml` has never passed, and it has to be fixed before the release gate can work.
- Anyone who can reach a preview's port can open it, which is why previews only ever hold synthetic data.
