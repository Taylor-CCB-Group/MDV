import "microtip/microtip.css";
import "../src/css/fontawesome-5.15.3/all.css";
import 'nouislider/dist/nouislider.min.css'
import "../src/utilities/css/ContextMenu.css";
import "../src/charts/css/charts.css";
import "../src/webgl/css/wgl2di.css";
import "../src/table/css/slickgrid.css";

import ChartManager from "../src/charts/ChartManager";


const dataSources=[{
    size:2735, //number of rows
    name:"cells", //name given here must be consistent with that in config and dataLoader
    //columns have:-
    //*datatype* - can be integer, double , text(<257 categories)
    //or unique (>256 categories e.g barcode)
    //*field* - this is how the column is referenced -in csv and tsv files it must
    //correspond to the headers or be the keys in a json file
    //*name* - the column label shown to the user - if no field parameter is present 
    //then the field will take on the same value as the name (as below)
    columns:[
    {
        "name":"cell_id",
        "datatype":"unique"
    },
    {
        "name":"UMAP_3d_1",
        "datatype":"double"
    },
    {
        "name":"UMAP_3d_2",
        "datatype":"double"
    },
    {
        "name":"UMAP_3d_3",
        "datatype":"double"
    },
    {
        "name":"gene_HCCAT5",
        "datatype":"double"
    },
    {
        "name":"gene_PART1",
        "datatype":"double"
    },
    {
        "name":"gene_THBD",
        "datatype":"double"
    },
    {
        "name":"gene_SIRPA",
        "datatype":"double"
    },
    {
        "name":"gene_TLR4",
        "datatype":"double"
    }, 				
    {
        "name":"Clusters",
        "datatype":"text"
    }],
    columnGroups:[
        {
            name:"genes",
            columns:["gene_HCCAT5","gene_PART1","gene_THBD","gene_SIRPA","gene_TLR4"]
        }
    ]
}];

//simple dataLoader - just specifies a single csv file
const dataLoader = {
    files:[
        {
            url:"data/pbmc3k.tsv", //url of file
            dataSource:"cells", //datasource (as specified above)
            type:"tsv" //can be tsv,csv or json
        }
    ]
};


 //config - single view
const config={
    only_view:{
        //the initial charts to load
        initialCharts:{
            "cells":[
            {
                type:"wgl_3d_scatter_plot",
                param:["UMAP_3d_1","UMAP_3d_2","UMAP_3d_3"],
                color_by:"Clusters",
                size:[400,600],
                position:[10,10]
            },
            {
                type:"row_chart",
                param:"Clusters",
                id:"my-clusters",
                size:[200,390],
                position:[420,220],
            }
        
            ]
        }
    }
};
//create the app with the above parameters
const cm = new ChartManager("app1",dataSources,dataLoader,config);