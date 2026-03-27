import click
import scanpy as sc
import mudata as mu
import shutil
import zipfile
import os
import json
from os.path import exists, join, basename
from .conversions import (
    convert_scanpy_to_mdv,
    convert_mudata_to_mdv,
    convert_vcf_to_mdv,
    merge_projects,
)
from .spatial.conversion import (
    convert_spatialdata_to_mdv,
    SpatialDataConversionArgs,
)

def zip_and_remove(folder):
    """Zip a directory and delete the original."""
    if not exists(folder):
        raise FileNotFoundError(f"{folder} not found")

    zip_filename = f"{folder.rstrip(os.sep)}.zip"
    click.echo(f"Zipping folder: {folder} → {zip_filename}")

    with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(folder):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, start=folder)
                zipf.write(file_path, arcname)

    click.echo(f"Zip created: {zip_filename}")
    click.echo(f"Removing original folder: {folder}")
    shutil.rmtree(folder)
    click.echo("Done.")

@click.group()
def cli():
    """MDV Tools CLI"""
    pass

@cli.command()
@click.argument('folder')
@click.argument('scanpy_object')
@click.option('--max_dims', default=3, help='Maximum number of dimensions to include from dimensionality reductions.')
@click.option('--delete_existing', is_flag=True, help='Delete existing project data.')
@click.option('--label', default='', help='Prefix to add to datasource names and metadata columns.')
@click.option('--chunk_data', is_flag=True, help='Transpose and flatten in chunks to save memory.')
@click.option('--add_layer_data', is_flag=True, default=True, help='Add layer data (log values etc.).')
@click.option('--gene_identifier_column', default=None, help='Gene column for identification.')
@click.option('--zip', 'zip_output', is_flag=True, help='Zip the output folder and delete the original.')
@click.option('--chatmdv', is_flag=True, help='Include the original Scanpy .h5ad file in the zipped project.')
def convert_scanpy(folder, scanpy_object, max_dims, delete_existing, label, chunk_data, add_layer_data, gene_identifier_column, zip_output, chatmdv):
    """Convert Scanpy AnnData object to MDV format."""
    adata = sc.read_h5ad(scanpy_object)
    convert_scanpy_to_mdv(folder, adata, max_dims, delete_existing, label, chunk_data, add_layer_data, gene_identifier_column)

    if chatmdv:
        dest_path = os.path.join(folder, basename(scanpy_object))
        shutil.copy(scanpy_object, dest_path)
        click.echo(f"Included original scanpy file: {dest_path}")

    if zip_output:
        zip_and_remove(folder)

@cli.command()
@click.argument('folder')
@click.argument('mudata_object')
@click.option('--max_dims', default=3, help='Maximum number of dimensions to include from dimensionality reductions.')
@click.option('--delete_existing', is_flag=True, help='Delete existing project data.')
@click.option('--chunk_data', is_flag=True, help='Transpose and flatten in chunks to save memory.')
@click.option('--zip', 'zip_output', is_flag=True, help='Zip the output folder and delete the original.')
def convert_mudata(folder, mudata_object, max_dims, delete_existing, chunk_data, zip_output):
    """Convert MuData object to MDV format."""
    mdata = mu.read(mudata_object)
    convert_mudata_to_mdv(folder, mdata, max_dims, delete_existing, chunk_data)
    if zip_output:
        zip_and_remove(folder)

@cli.command()
@click.argument('folder')
@click.argument('vcf_filename')
@click.option('--zip', 'zip_output', is_flag=True, help='Zip the output folder and delete the original.')
def convert_vcf(folder, vcf_filename, zip_output):
    """Convert VCF file to MDV format."""
    convert_vcf_to_mdv(folder, vcf_filename)
    if zip_output:
        zip_and_remove(folder)

@cli.command()
@click.argument('folder')
@click.option('--port', default=5050, help='Port to serve on.')
def serve(folder, port):
    """Serve MDV project."""
    from .serverlite import serve_project
    from .mdvproject import MDVProject
    if not exists(folder):
        raise FileNotFoundError(f"{folder} not found")
    ds_path = join(folder, "datasources.json")
    if not exists(ds_path):
        raise FileNotFoundError(f"{folder} does not contain a valid MDV project.")
    serve_project(MDVProject(folder), port=port)


@cli.command("merge-project")
@click.argument("base_project")
@click.argument("extra_project")
@click.option(
    "--prefix",
    default=None,
    help="Prefix applied to imported datasource names. "
    "If omitted, a prefix derived from EXTRA_PROJECT is used.",
)
@click.option(
    "--view-prefix",
    default=None,
    help="Prefix applied to imported view names. Defaults to the datasource prefix.",
)
def merge_project(base_project, extra_project, prefix, view_prefix):
    """Merge an existing MDV project into another project."""
    merge_projects(base_project, extra_project, prefix=prefix, view_prefix=view_prefix)
    click.echo(f"Merged '{extra_project}' into '{base_project}'.")


def _parse_csv_values(raw: str) -> list[str]:
    values = [item.strip() for item in raw.split(",")]
    return [value for value in values if value]


def _parse_rename_map(raw: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    if not raw.strip():
        return mapping
    for pair in raw.split(","):
        entry = pair.strip()
        if not entry:
            continue
        if ":" not in entry:
            raise click.BadParameter(
                f"Invalid mapping '{entry}'. Expected format old:new."
            )
        old, new = entry.split(":", 1)
        old = old.strip()
        new = new.strip()
        if not old or not new:
            raise click.BadParameter(
                f"Invalid mapping '{entry}'. Both old and new names are required."
            )
        mapping[old] = new
    return mapping


def _print_impact_summary(report: dict) -> None:
    impact = report.get("impactReport") or {}
    charts = impact.get("charts") or []
    state_refs = impact.get("stateRefs") or []
    ds_refs = impact.get("datasourceRefs") or []

    if not charts and not state_refs and not ds_refs:
        click.echo("Impact summary: no affected references found.")
        return

    click.echo("Impact summary:")
    if charts:
        by_view: dict[str, int] = {}
        by_type: dict[str, int] = {}
        for item in charts:
            view = item.get("view") or "<unknown>"
            by_view[view] = by_view.get(view, 0) + 1
            ctype = (item.get("chart") or {}).get("type") or "<unknown>"
            by_type[ctype] = by_type.get(ctype, 0) + 1
        click.echo(f"- charts affected: {len(charts)}")
        click.echo(f"- views affected: {len(by_view)}")
        top_views = sorted(by_view.items(), key=lambda x: (-x[1], x[0]))[:5]
        click.echo(f"- top views: {', '.join([f'{k}({v})' for k,v in top_views])}")
        top_types = sorted(by_type.items(), key=lambda x: (-x[1], x[0]))[:5]
        click.echo(f"- top chart types: {', '.join([f'{k}({v})' for k,v in top_types])}")
    if state_refs:
        click.echo(f"- state refs affected: {len(state_refs)}")
    if ds_refs:
        click.echo(f"- datasource nested refs affected: {len(ds_refs)}")

@cli.command("list-columns")
@click.argument("folder")
@click.option("--datasource", required=True, help="Datasource name to inspect.")
@click.option("--not-used", "not_used", is_flag=True, help="Return only columns not referenced in project configs.")
@click.option("--used-only", "used_only", is_flag=True, help="Return only columns referenced in project configs.")
@click.option("--json", "as_json", is_flag=True, help="Print full JSON usage report.")
@click.option("--one-per-line", is_flag=True, help="Print selected fields one per line.")
def list_columns(folder, datasource, not_used, used_only, as_json, one_per_line):
    """List columns in a datasource, optionally filtering to not-used fields."""
    from .mdvproject import MDVProject

    if not_used and used_only:
        raise click.BadParameter("--not-used and --used-only are mutually exclusive.")
    project = MDVProject(folder)
    report = project.list_columns(datasource, not_used=not_used)
    if used_only:
        fields = report["usedFields"]
    else:
        fields = report["fields"]
    if as_json:
        click.echo(json.dumps(report, indent=2, sort_keys=True))
        return
    if one_per_line:
        click.echo("\n".join(fields))
        return
    click.echo(",".join(fields))


@cli.command("drop-columns")
@click.argument("folder")
@click.option("--datasource", required=False, help="Datasource name to update.")
@click.option("--fields", required=False, help="Comma-separated field ids to drop.")
@click.option("--strict/--no-strict", default=True, show_default=True)
@click.option("--dry-run", is_flag=True, help="Preview changes without writing.")
@click.option("--cleanup", is_flag=True, help="Rewrite project JSON to remove references to dropped fields.")
@click.option("--backup", is_flag=True, help="Backup JSON files before rewriting.")
@click.option("--impact-summary/--no-impact-summary", default=True, show_default=True, help="Print impact summary on dry-run.")
@click.option("--restore-backup", "restore_backup", is_flag=True, help="Restore latest JSON backup and exit.")
@click.option("--restore-backup-timestamp", default=None, help="Restore a specific JSON backup timestamp and exit.")
@click.option("--continue-after-restore", is_flag=True, help="After restoring a backup, continue with the requested operation.")
@click.option(
    "--hard",
    is_flag=True,
    help="Perform physical HDF5 dataset deletion in addition to metadata updates.",
)
def drop_columns(folder, datasource, fields, strict, dry_run, cleanup, backup, impact_summary, restore_backup, restore_backup_timestamp, continue_after_restore, hard):
    """Drop multiple columns from a datasource."""
    from .mdvproject import MDVProject

    if hard:
        click.echo("WARNING: --hard will mutate HDF5 datasets.")
    if backup and not cleanup:
        raise click.BadParameter("--backup requires --cleanup.")
    project = MDVProject(folder)
    if restore_backup or restore_backup_timestamp:
        report = project.restore_json_backup(
            timestamp=restore_backup_timestamp,
            dry_run=dry_run,
        )
        click.echo(json.dumps(report, indent=2, sort_keys=True))
        if not continue_after_restore:
            return
    if not datasource:
        raise click.BadParameter("Missing option '--datasource'.")
    if not fields:
        raise click.BadParameter("Missing option '--fields'.")
    parsed_fields = _parse_csv_values(fields)
    if not parsed_fields:
        raise click.BadParameter("At least one field must be provided.", param_hint="fields")
    try:
        report = project.drop_columns(
            datasource,
            parsed_fields,
            strict=strict,
            hard=hard,
            dry_run=dry_run,
            cleanup=cleanup,
            backup=backup,
        )
        if dry_run and impact_summary:
            _print_impact_summary(report)
        click.echo(json.dumps(report, indent=2, sort_keys=True))
    except ValueError as e:
        details = getattr(e, "details", None)
        if details:
            click.echo(json.dumps(details, indent=2, sort_keys=True))
        raise


@cli.command("rename-columns")
@click.argument("folder")
@click.option("--datasource", required=False, help="Datasource name to update.")
@click.option(
    "--map",
    "rename_map",
    required=False,
    help="Comma-separated rename map in old:new format.",
)
@click.option(
    "--field-id",
    is_flag=True,
    help="Rename internal field ids instead of display labels.",
)
@click.option("--strict/--no-strict", default=True, show_default=True)
@click.option("--dry-run", is_flag=True, help="Preview changes without writing.")
@click.option("--cleanup", is_flag=True, help="Rewrite project JSON to apply field-id renames consistently.")
@click.option("--backup", is_flag=True, help="Backup JSON files before rewriting.")
@click.option("--impact-summary/--no-impact-summary", default=True, show_default=True, help="Print impact summary on dry-run.")
@click.option("--restore-backup", "restore_backup", is_flag=True, help="Restore latest JSON backup and exit.")
@click.option("--restore-backup-timestamp", default=None, help="Restore a specific JSON backup timestamp and exit.")
@click.option("--continue-after-restore", is_flag=True, help="After restoring a backup, continue with the requested operation.")
@click.option(
    "--hard",
    is_flag=True,
    help="Required for operations that mutate HDF5 datasets.",
)
def rename_columns(folder, datasource, rename_map, field_id, strict, dry_run, cleanup, backup, impact_summary, restore_backup, restore_backup_timestamp, continue_after_restore, hard):
    """Rename multiple columns on a datasource."""
    from .mdvproject import MDVProject

    project = MDVProject(folder)
    if restore_backup or restore_backup_timestamp:
        report = project.restore_json_backup(
            timestamp=restore_backup_timestamp,
            dry_run=dry_run,
        )
        click.echo(json.dumps(report, indent=2, sort_keys=True))
        if not continue_after_restore:
            return
    if not datasource:
        raise click.BadParameter("Missing option '--datasource'.")
    if not rename_map:
        raise click.BadParameter("Missing option '--map'.")
    renames = _parse_rename_map(rename_map)
    if field_id and not hard:
        raise click.BadParameter("--field-id requires --hard.")
    if backup and not cleanup:
        raise click.BadParameter("--backup requires --cleanup.")
    if hard:
        click.echo("WARNING: --hard will mutate HDF5 datasets.")
    report = project.rename_columns(
        datasource,
        renames,
        rename_field=field_id,
        strict=strict,
        hard=hard,
        dry_run=dry_run,
        cleanup=cleanup,
        backup=backup,
    )
    if dry_run and impact_summary:
        _print_impact_summary(report)
    click.echo(json.dumps(report, indent=2, sort_keys=True))

@cli.command("convert-spatial")
@click.argument('output_folder')
@click.argument('spatialdata_path')
@click.option('--preserve-existing', 'preserve_existing', is_flag=True, help='Preserve existing project data.')
@click.option('--link', is_flag=True, help='Symlink to the original SpatialData objects.')
@click.option('--output_geojson', is_flag=True, help='Output geojson for each region (this feature to be deprecated in favour of spatialdata.js layers with shapes).')
@click.option('--serve', is_flag=True, help='Serve the project after conversion.')
def convert_spatial(spatialdata_path, output_folder, preserve_existing, link, output_geojson, serve):
    """Convert SpatialData objects to MDV format."""
    import tempfile
    with tempfile.TemporaryDirectory() as temp_folder:
        args = SpatialDataConversionArgs(
            spatialdata_path=spatialdata_path,
            output_folder=output_folder,
            preserve_existing=preserve_existing,
            temp_folder=temp_folder,
            link=link,
            output_geojson=output_geojson,
            serve=serve,
        )
        convert_spatialdata_to_mdv(args)

if __name__ == '__main__':
    cli()
