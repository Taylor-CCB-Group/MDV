import click
from os.path import exists,join
from .conversions import convert_scanpy_to_mdv, convert_mudata_to_mdv, convert_vcf_to_mdv

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
def convert_scanpy(folder, scanpy_object, max_dims, delete_existing, label, chunk_data, add_layer_data, gene_identifier_column):
    """Convert Scanpy AnnData object to MDV format."""
    convert_scanpy_to_mdv(folder, scanpy_object, max_dims, delete_existing, label, chunk_data, add_layer_data, gene_identifier_column)

@cli.command()
@click.argument('folder')
@click.argument('mudata_object')
@click.option('--max_dims', default=3, help='Maximum number of dimensions to include from dimensionality reductions.')
@click.option('--delete_existing', is_flag=True, help='Delete existing project data.')
@click.option('--chunk_data', is_flag=True, help='Transpose and flatten in chunks to save memory.')
def convert_mudata(folder, mudata_object, max_dims, delete_existing, chunk_data):
    """Convert MuData object to MDV format."""
    convert_mudata_to_mdv(folder, mudata_object, max_dims, delete_existing, chunk_data)

@cli.command()
@click.argument('folder')
@click.argument('vcf_filename')
def convert_vcf(folder, vcf_filename):
    """Convert VCF file to MDV format."""
    convert_vcf_to_mdv(folder, vcf_filename)



@cli.command()
@click.argument('folder')
def serve(folder):
    """Serve MDV project."""
    from .serverlite import serve_project
    from .mdvproject import MDVProject
    if not exists(folder):
            raise FileNotFoundError(f"{folder} not found")
    #could do more extensive check
    ds_path = join(folder, "datasources.json")
    if not exists(ds_path):
        raise FileNotFoundError(f"{folder} does not contain a valid MDV project.")
    serve_project(MDVProject(folder))

if __name__ == '__main__':
    cli()