import scanpy
import mudata
from .mdvproject import MDVProject

def convert_scanpy_to_mdv(folder,scanpy_object,max_dims=3):
    mdv =  MDVProject(folder)

    #create datasources 'cells'
    cell_table = scanpy_object.obs
    cell_table["cell_id"] = cell_table.index
    #add any dimension reduction
    _add_dims(cell_table,scanpy_object.obsm,max_dims)
    mdv.add_datasource("cells",cell_table)

     #create two datasources 'genes'
    gene_table = scanpy_object.var
    gene_table["name"]=gene_table.index
    _add_dims(gene_table,scanpy_object.varm,max_dims)
    mdv.add_datasource("genes",gene_table)

    #link the two datasets
    mdv.add_rows_as_columns_link("cells","genes","name","Gene Expr")

    #add the gene expression 
    mdv.add_rows_as_columns_subgroup("cells","genes","gs",scanpy_object.X,name="gene_scores",label="Gene Scores")

    #create a default view
    mdv.set_view("default",{
        "initialCharts":{
            "cells":[],
            "genes":[]
        }
    },True)
    return mdv

    
def convert_mudata_to_mdv(folder,mudata_object,max_dims=3):
    md=mudata_object
    p= MDVProject(folder)
    #add the cells
    p.add_datasource("cells",md.obs)

    for mod in md.mod.keys():   
        mdata = md.mod[mod]
        if hasattr(mdata.X,"value"):
            p.add_datasource(mod,mdata.var)
            p.add_column(mod,{"name":"name","datatype":"unique"},mdata.var.index)
            p.add_rows_as_columns_link("cells",mod,"name",mod)
            matrix = mdata.X.value.transpose().tocsr()
            p.add_rows_as_columns_subgroup("cells",mod,mod+"_expr",matrix,sparse=True)



    
def _add_dims(table,dims,max_dims):
    for dname in  dims:
        dm= dims[dname].transpose()
        md= min(max_dims,dm.shape[0])
        for n in range(md):
            table[f"{dname}_{n+1}"]=dm[n]







