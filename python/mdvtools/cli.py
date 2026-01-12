import click
import scanpy as sc
import mudata as mu
import shutil
import zipfile
import os
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

@cli.command("convert-spatial")
@click.argument('spatialdata_path')
@click.argument('output_folder')
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
