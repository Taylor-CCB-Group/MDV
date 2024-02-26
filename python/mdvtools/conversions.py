import scanpy
import mudata
import scipy
from os.path import join,split
import os
from .mdvproject import MDVProject
import pandas as pd
import numpy as np
import json
import gzip

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
            p.set_column(mod,{"name":"name","datatype":"unique"},mdata.var.index)
            p.add_rows_as_columns_link("cells",mod,"name",mod)
            matrix = scipy.sparse.csr_matrix(matrix).transpose().tocsr()
            p.add_rows_as_columns_subgroup("cells",mod,mod+"_expr",matrix,sparse=True)

    return p


def create_regulamentary_project(folder,output,name,
                bigwig_folder="",
                bed_folder="",
                openchrome="DNase",
                openchrome_color="#eb9234",
                marks= ["H3K4me1","H3K4me3","H3K27ac","CTCF"],
                mark_colors=["#349beb","#3aeb34","#c4c41f","#ab321a"],
                genome = "hg38"):
   

    #get the template dir
    tdir = join(split(os.path.abspath(__file__))[0],"templates")
    p = MDVProject(output)
    mdvfile  =  join(folder,"08_REgulamentary","mlv_REgulamentary.csv")
    mdv = pd.read_csv(mdvfile,sep="\t")
    p.add_datasource("elements",mdv)
  
   
    #add the tracks
    all_names=  [openchrome]+marks
    all_colors= [openchrome_color]+mark_colors
    if bigwig_folder.startswith("http"):
        bigwigs = [f"{bigwig_folder}/{name}_{x}_{genome}.bw" for x in all_names]
        beds= [f"{bed_folder}/{name}_{x}_{genome}.bb" for x in all_names]
    else:
        pass
    #get the reference
  
    default_tracks=[]
    for n in range(len(all_names)):
        default_tracks.append({
            "url":bigwigs[n],
            "short_label":all_names[n]+" cov",
            "color": all_colors[n],
            "height":60,
            "track_id":"coverage_"+all_names[n]
        })
        default_tracks.append({
            "url":beds[n],
            "short_label":all_names[n] + " peaks",
            "color": all_colors[n],
            "height":15,
            "featureHeight":5,
            "track_id":"peaks_"+all_names[n]
        })
       
       
    extra_params={
        "default_tracks":default_tracks,
        "default_parameters":{
            "view_margins": {
                    "type": "fixed_length",
                    "value": 4000
            },  
            "color_by": "RE",
            "color_legend": {
                "display": False,
                "pos": [
                5,
                5
                ]
            },
            "feature_label": "RE",
            "view_margins": {
                "type": "fixed_length",
                "value": 4000
            }
        }
    }
    p.add_genome_browser("elements",["chromosome","start","end"], 
                         extra_params=extra_params)
    p.add_refseq_track("elements",genome)
    
    _create_dt_heatmap(folder,p,mdv,marks)
    
    
    view = json.load(open(join(tdir,"views","regulamentary.json")))
    gb = p.get_genome_browser("elements")
    gb.update({
        "id":"browser",
        "size":[389,550],
        "position":[1046,8],
        "title":"elements"

    })
    view.append(gb)

    p.set_view("default",{
        "initialCharts":{
            "elements":view,
        }
    })
    return p

    

def _create_dt_heatmap(folder,project,mdv,marks):
    #get the order of the regions in the matrix
    mdv_i = mdv.set_index(["chromosome","start","end"])
    order = pd.read_csv(join(folder,"04_sort_regions","sort_union.bed"),sep="\t",header=None)
    order= order.set_index([0,1,2])
    order= order[order.index.isin(mdv_i.index)]
    mdv_i=mdv.reset_index()
    mdv_i['row_position'] = mdv_i.index
    arr = order.index.map(mdv_i.set_index(['chromosome', 'start', 'end'])['row_position'])
    arr = np.array(arr,dtype=np.int32)

    #add the binary data
    output = join(project.dir,"binarydata","elements")
    os.makedirs(output,exist_ok=True)
    hm = gzip.open(join(output,"heat_map.gz"),"wb")
    hm.write(arr.tobytes())

    #get the matrix and drop last two columns
    mt = pd.read_csv(join(folder,"09_metaplot","matrix.csv"),sep="\t")
    mt = mt.iloc[:,:1600]
    #flatten and halve in size (mean of adjacent values)
    arr = mt.values.flatten()
    arr  = arr.reshape(-1,2)
    arr= np.mean(arr,axis=1)
    #normalize and clip to 255 (1 byte per value)
    max = np.percentile(arr, 99.99)
    arr = ((arr/max)*255)
    arr= np.clip(arr,0,255)
    arr=arr.astype(np.uint8)
    hm.write(arr.tobytes())
    hm.close()
    #add the metdata
    md = project.get_datasource_metadata("elements")
    md["deeptools"]={
       "maps":{
         "default" :{
            "data":"heat_map",
             "rows":md["size"],
             "cols":800,
             "groups":marks,
            "max_scale":max
         }
      }
    }
    md["binary_data_loader"]=True
    project.set_datasource_metadata(md)
       
     


   

    
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







