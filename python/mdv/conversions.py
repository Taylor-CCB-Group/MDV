import scanpy
import mudata
import scipy

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
    _add_dims(md.obs,md.obsm,max_dims)
    p.add_datasource("cells",md.obs)

    for mod in md.mod.keys():   
        mdata = md.mod[mod]
        matrix = mdata.X.value if hasattr(mdata.X,"value") else mdata.X
        if matrix.shape[1] !=0:
            p.add_datasource(mod,mdata.var)
            p.add_column(mod,{"name":"name","datatype":"unique"},mdata.var.index)
            p.add_rows_as_columns_link("cells",mod,"name",mod)
            matrix = scipy.sparse.csr_matrix(matrix).transpose().tocsr()
            p.add_rows_as_columns_subgroup("cells",mod,mod+"_expr",matrix,sparse=True)

    return p

    
def _add_dims(table,dims,max_dims):
    if len(dims.keys())==0:
        return
    for dname in  dims:
        if len(dims[dname].shape)==1:
            continue
        dm= dims[dname].transpose()
        md= min(max_dims,dm.shape[0])
        for n in range(md):
            table[f"{dname}_{n+1}"]=dm[n]







