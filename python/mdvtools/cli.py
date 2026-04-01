import click
import shutil
import zipfile
import os
from os.path import exists, join, basename

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
@click.option('--obs-datasource-name', default='cells', show_default=True, help='Datasource name for observations.')
@click.option('--var-datasource-name', default='genes', show_default=True, help='Datasource name for variables.')
@click.option('--chunk_data', is_flag=True, help='Transpose and flatten in chunks to save memory.')
@click.option('--add_layer_data', is_flag=True, default=True, help='Add layer data (log values etc.).')
@click.option('--gene_identifier_column', default=None, help='Gene column for identification.')
@click.option('--link-name-column', default=None, help='Variable datasource column used as the rows_as_columns name_column.')
@click.option('--compute-x-umap', 'compute_x_umap', is_flag=True, help='Compute neighbors, UMAP, and Leiden clusters directly from adata.X before export.')
@click.option('--leiden-resolution', default=1.0, type=float, show_default=True, help='Leiden resolution used with --compute-x-umap.')
@click.option('--zip', 'zip_output', is_flag=True, help='Zip the output folder and delete the original.')
@click.option('--chatmdv', is_flag=True, help='Include the original Scanpy .h5ad file in the zipped project.')
def convert_scanpy(folder, scanpy_object, max_dims, delete_existing, label, obs_datasource_name, var_datasource_name, chunk_data, add_layer_data, gene_identifier_column, link_name_column, compute_x_umap, leiden_resolution, zip_output, chatmdv):
    """Convert Scanpy AnnData object to MDV format."""
    import scanpy as sc
    from .conversions import convert_scanpy_to_mdv

    adata = sc.read_h5ad(scanpy_object)
    convert_scanpy_to_mdv(
        folder,
        adata,
        max_dims,
        delete_existing,
        label,
        obs_datasource_name,
        var_datasource_name,
        chunk_data,
        add_layer_data,
        gene_identifier_column,
        link_name_column=link_name_column,
        compute_x_umap=compute_x_umap,
        leiden_resolution=leiden_resolution,
    )

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
    import mudata as mu
    from .conversions import convert_mudata_to_mdv

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
    from .conversions import convert_vcf_to_mdv

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
    from .conversions import merge_projects

    merge_projects(base_project, extra_project, prefix=prefix, view_prefix=view_prefix)
    click.echo(f"Merged '{extra_project}' into '{base_project}'.")

@cli.command("convert-spatial")
@click.argument('output_folder')
@click.argument('spatialdata_path')
@click.option('--batch', is_flag=True, help='Convert all child SpatialData stores in the given directory.')
@click.option('--preserve-existing', 'preserve_existing', is_flag=True, help='Preserve existing project data in the output folder instead of recreating it.')
@click.option('--link', is_flag=True, help='Symlink the original SpatialData source into the project instead of copying it.')
@click.option('--output_geojson/--no-output_geojson', 'output_geojson', default=True, help='Write transformed GeoJSON region files into the project images directory.')
@click.option('--density', is_flag=True, help='Include density fields for gene expression in the default spatial view.')
@click.option('--serve', is_flag=True, help='Serve the generated project after conversion.')
@click.option('--obs-datasource-name', default='cells', show_default=True, help='Datasource name for observations.')
@click.option('--var-datasource-name', default='genes', show_default=True, help='Datasource name for variables.')
@click.option('--link-name-column', default=None, help='Variable datasource column used as the rows_as_columns name_column.')
@click.option('--compute-x-umap', 'compute_x_umap', is_flag=True, help='Compute neighbors, UMAP, and Leiden clusters from merged adata.X before export.')
@click.option('--leiden-resolution', default=1.0, type=float, show_default=True, help='Leiden resolution used with --compute-x-umap.')
@click.option('--verbose', is_flag=True, help='Show detailed per-dataset conversion output, transform decisions, and merged summaries.')
def convert_spatial(spatialdata_path, output_folder, batch, preserve_existing, link, output_geojson, density, serve, obs_datasource_name, var_datasource_name, link_name_column, compute_x_umap, leiden_resolution, verbose):
    """Convert one SpatialData store to MDV format, or use --batch for a directory of stores."""
    import tempfile
    from .spatial.conversion import convert_spatialdata_to_mdv, SpatialDataConversionArgs

    with tempfile.TemporaryDirectory() as temp_folder:
        args = SpatialDataConversionArgs(
            spatialdata_path=spatialdata_path,
            output_folder=output_folder,
            batch=batch,
            preserve_existing=preserve_existing,
            temp_folder=temp_folder,
            link=link,
            output_geojson=output_geojson,
            density=density,
            serve=serve,
            verbose=verbose,
            obs_datasource_name=obs_datasource_name,
            var_datasource_name=var_datasource_name,
            link_name_column=link_name_column,
            compute_x_umap=compute_x_umap,
            leiden_resolution=leiden_resolution,
        )
        convert_spatialdata_to_mdv(args)

if __name__ == '__main__':
    cli()
