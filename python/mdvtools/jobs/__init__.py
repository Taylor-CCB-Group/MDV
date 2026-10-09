# Name of the per-project jobs directory holding the durable owner-side records.
# It lives *inside* the project (`<project>/jobs/records/`) so the catalog's top-level
# scan (mdv_server_app.serve_projects_from_filesystem / mdv_desktop) never mistakes it
# for a project, and so provenance travels with the project. The exporter skips it
# (mdvproject.convert_to_static_page). The ephemeral per-job scratch is a *separate*
# root outside the project (ADR-0007); the two are linked only by job_id.
JOBS_DIRNAME = "jobs"
