import ChartManager from "../charts/ChartManager";
import pako from "pako";

async function scLoadMetaData(folder){
    const proj={
        data:{},
        permission:"view"
    }
    let resp = await fetch(`${folder}/datasources.json`);
    proj.data.dataSources =  await resp.json();
    resp = await fetch(`${folder}/state.json`);
    proj.data.state =  await resp.json();
    resp = await fetch(`${folder}/hyperion_config.json`);
    proj.data.hyperion_config =  await resp.json();
    scLoadData(proj,folder);

}


function scLoadData(proj,local){
    let config = proj.data.state;
    config["permission"]=proj.permission;
    //legacy
    if (config.initialCharts){
        config.only_view={
            initialCharts:config.initialCharts
        }
        delete config.initialCharts;
    }
    if (!proj.data.hyperion_config.default_images){
        proj.data.hyperion_config.default_images={};
    }

    const urlParams = new URLSearchParams(window.location.search);
    const view = urlParams.get("view");
    if (config.all_views && view && config.all_views.indexOf(view) !==-1){
        config["initial_view"]=view;
    }
    const dataLoader= 
        local?
        {
            function:scGetLocalDataLoader(local,proj.data.dataSources),
            viewLoader:scGetLocalViewLoader(local)

        }:
        {
            function:scDataLoader,
            viewLoader:scViewLoader
        };

    const listener = (type,cm,data)=>{
        scListener(type,cm,data,proj);
    }
    const cm = new ChartManager("app-holder",proj.data.dataSources,dataLoader,config,listener); 
    cm.addMenuIcon("_main","fas fa-cogs","re-cluster selected cells",()=>scRecluster(cm,proj.data.hyperion_config));
    cm._ssfilter=_cm.dsIndex["cells"].dataStore.getDimension("category_dimension");

    const urlbase= local?`${local}/images`:`/data/${proj.genome}/projects/${proj.id}/images`;
    cm._urlbase=urlbase;
    cm.hyperion_config=proj.data.hyperion_config;
    for (let  name in cm.dsIndex){
        const ds=  cm.dsIndex[name];
        if (ds.custom.spatial_stats){
            scAddSSListener(ds,urlbase,cm,proj.data.hyperion_config);
            ds._ssfilter=ds.dataStore.getDimension("category_dimension");
        }
    }

}

function scAddSSListener(ds,urlbase,cm,hyp_conf){
    ds.dataStore.addListener("ss_listen",(type)=>{
        if (type==="data_highlighted"){
            scShowCellAdjacency(ds,urlbase,cm,hyp_conf);
        }
    });
}

function scRecluster(cm,hyp_conf){
    const ds = cm.getDataSource("cells");
    const filters = cm.getAllFilters("cells");
    let f_text= "The data will be reclustered using only the currently filtered cells and the markers selected below<br>";
    f_text+=filters.length===0?"<b>No Filters</b><br>":"<b>Current Filters<b>:<br>";

    const markers = hyp_conf.all_markers || ds.getColumnGroup("Markers");
    const cbm = markers.map(x=>{
        return {
            id:x,
            label:x,
            type:"checkbox",
            current_value:true
        }
    })

    const r_choices=[
        ["All","all"],

    ]

    if (hyp_conf.clustering){
        for (c in hyp_conf.clustering){
            r_choices.push([c,c])
        }
    }


    for (let f of filters){
        for (let i in f){
        const col = ds.columnIndex[i];
        if (col.datatype!="text"){
            f[i]=f[i].map(x=>x.toPrecision(4));
        }
        f_text+=`${col.name}:${f[i].join(col.datatype==="text"?",":"-")}<br>`;
        }
    }
    f_text+="<b>Number of Cells<b>:"+ds.filterSize;
    cm.showCustomDialog({
    text:f_text,
    title:"Re-Cluster",
    maxHeight:400,
    buttons:[{
        text:"submit",
        method:(c)=>{
            const mks= [];
            for (let m in c){
                if (markers.indexOf(m)!==-1){
                    if (c[m]){
                        mks.push(m);
                    }
                }
            }
            scSubmitClustering(mks,c["name"],ds,cm);
        }
        
        }
        
    ],
    controls:[
        {
            type:"text",
            label:"Name",
            description:"A short memorable label that will used as a prefix for views and columns",
            id:"name"
        },
        {
            type:"radiobuttons",
            choices:r_choices,
            label:"Markers",
            id:"mcs",
            func:(v,dialog)=>{
                for (let m of markers){
                    const cb = dialog.controls[m];
                    val = true;
                    if (v!=="all"){
                        val = hyp_conf.clustering[v].markers.indexOf(m) !== -1
                    }
                    cb.checked=val;
                    dialog.controlValues[m]=val;
                }
            }
        }		
        ].concat(cbm)
    })
}

function scSubmitClustering(markers,name,dataStore,cm){
    cm._getColumnsThen("cells",["sample_id"],()=>{
        const f= dataStore.getFilter();
        const ids= [];
        const all_samples  = dataStore.columnIndex["sample_id"].data;
        const samples= new Set()
        for (let i=0;i<f.length;i++){
            if (f[i]===0){
                samples.add(all_samples[i]);
                ids.push(i)
            }	
        }
        executeProjectAction("recluster",{
            ids:ids,
            name:name,
            markers:markers,
            filters:cm.getAllFilters("cells"),
            samples:Array.from(samples)
            }).then(resp=>{
                const msg = resp.success?
                    `Your Clustering job has been queued - You will receive an email when it is complete`:
                    `There was a problem, please contact an administrator`;
                
                cm.showCustomDialog({
                    title:"Message",
                    text:msg
                    })
            });
    });
}

function scAddImage(sample,urlbase,cm,hconf,iconf={}){
try{
    const info =hconf. images[sample];
    const urls= info.map(x=>{
        return [x[0],`${urlbase}/${x[1]}.png`,x[2],x[3],x[4],x[5]];
    });
    
    let di = info[0];
    if (hconf.default_image){
        const ndi = info.find(x=>x[0]===hconf.default_image);
        if (ndi){
            di=ndi;
        }
    }
    const densParam = iconf.color_by || hconf.default_color_field;

    const chart = {
        type:"density_scatter_plot",
        title:`${sample}-${di[0]}`,
        background_filter:{
            column:"sample_id",
            category:sample
        },
        color_legend:{
            display:false
        },
        legend:`The backround image can be changed in the settings dialog (cog icon) .Also, the image can
        be resized and moved by pressing ctrl amd zooiming or panning with the mouse`,
        radius:3,
        param:hconf.position_fields.concat([densParam]),
        roi:{
            max_x:hconf.roi_sizes[sample][0],
            max_y:hconf.roi_sizes[sample][1]
            
        },
        background_image:{
            url:`${urlbase}/${di[1]}.png`,
            position:[di[4],di[5]],
            width:di[2],
            height:di[3],
            name:di[0]
        },
        image_choices:urls
    }
    if (iconf.color_by){
        chart.color_by= iconf.color_by;
    }
    else if (hconf.default_color_field){
        chart.color_by = hconf.default_color_field;
    }
    if (iconf.id){
        chart.id=iconf.id;
    }
    if (iconf.size){
        chart.size=iconf.size;
    }
    const d= iconf.details;
    if (d){
        chart.opacity=d.opacity;
        if (d.cats){
            chart.category1=d.cats[0];
            chart.category2= d.cats[1];
        }
        if (d.point_size){
            chart.radius=d.point_size;
        }
    }
    if (hconf.default_images && hconf.default_images[sample]){
        chart.background_image= hconf.default_images[sample];
    }
    cm.addChart("cells",chart,true);
}catch(err){
    console.log(err);
}
}


function scDataLoader(columns,dataSource,size){
    const toSend= {
    method:"get_column_data",
    args:{
        columns:columns,
        data_source:dataSource
    }
}
return new Promise((resolve,reject)=>{
  fetch("/meths/execute_project_action/"+project_id,
    {
        method:"POST",
        body:JSON.stringify(toSend),
        headers:{
            "Accept":"application/json,text/plain,*/*",
            "Content-Type":"application/json"
        }
    })

  .then(response=>response.arrayBuffer())
  .then(data=>{
    const dataList= [];
    let offset=0;
   
    for (let column of columns){
      let arrayType=Float32Array;
     let arr_len=size;
      let bytes=4;
      if (column.datatype=="unique" || column.datatype=="text"){
          arrayType=Uint8Array;
          bytes=1;
          if(column.datatype=="unique"){
            bytes=column.stringLength;
            arr_len=size*bytes;
          }
         
      }
      const len  = size*bytes; 
    let arr = null
    try{
          arr = new arrayType(data,offset,arr_len);
    }catch(e){
        console.log(column);
    }
      const sab = new SharedArrayBuffer(len);
      const new_arr =  new arrayType(sab)
      /*for (let n=0;n<len;n++){
          new_arr[n]=arr[n];
      }*/
      new_arr.set(arr,0);
      dataList.push({data:sab,field:column.field})
      offset+=len;
    }       
    resolve(dataList);
  }).catch(reject);   
});

}

function scShowCellAdjacency(ds,urlbase,cm,hconf){
const cds =ds.dataStore;
if (ds.custom.main_pivot!=="sample_id"){
    return;
}
const cat = ds.column_link_to.columns[1].link_to;
const hl = cds.getHighlightedData();
const settings = cm.viewData.spatial_stats[ds.name];
if (!hl){
    return;
}
    cm._getColumnsThen(ds.name,["sample_id","Cell Type 1","Cell Type 2"],()=>{
        const info = cds.getRowAsObject(hl[0],["sample_id","Cell Type 1","Cell Type 2"]);
        const cid = "ss_"+info["sample_id"];
        const adp = cm.getChart(cid);
        
        if (!adp){
            const details={
                opacity:settings.show_points?0.8:0,
                cats:settings.show_density?[info["Cell Type 1"],info["Cell Type 2"]]:null
            }
            scAddImage(info["sample_id"],urlbase,cm,hconf,{color_by:cat,id:cid,details:details});
        }
        cm._getColumnsThen("cells",[cat],()=>{
            if (adp ){
                if (settings.show_density){
                    adp.config.category1= info["Cell Type 1"];
                    adp.config.category2= info["Cell Type 2"];
                }
                adp.config.opacity=settings.show_points?0.8:0;
            }
            cm._ssfilter.filter("filterCategories",[cat],[info["Cell Type 1"],info["Cell Type 2"]])
        })
})
}

function scInitApp(cm,proj){

const hconf = proj.data.hyperion_config;
if (hconf){
    
    
    
        const chs = cm.getAllCharts("cells").filter(x=>x.config.background_image);
        for (let c of chs){
            scAddSetDefaultChoice(c,hconf);
        }
        
        const urlbase= `/data/${proj.genome}/projects/${proj.id}/images`;

        cm.addMenuIcon("cells","far fa-image","Add sample scatterplot with image",()=>{
            const choices=[];
            for (let i in hconf.images){
             choices.push({text:i,value:i})
        }
            cm.showCustomDialog({
                  title:"Add Image Scatter Plot",
                  text:"Select Sample:",
                  controls:[{
                    label:"",
                    id:"samples",
                    type:"dropdown",
                    items:choices
            
                  }],
                  buttons:[{
                    text:"Add",
                    method:c=>scAddImage(c.samples,urlbase,cm,hconf)
                  }]
            
            })
        });
        
        
        cm.addMenuIcon("cells","fas fa-file-upload","upload image",()=>{
            scShowUploadImageDialog(cm,hconf.images,urlbase)
        });

        cm.addMenuIcon("cells","far fa-images","bulk manipulate images",()=>{
            scShowBulkChangeImagesDialog(cm,hconf)
        });

        cm.addMenuIcon("cells","far fa-chart-bar","Add cell abundance chart",()=>{
            scShowCellAbundanceDialog(cm,hconf);

        });
        if (hconf.avivator){
            cm.addMenuIcon("cells","fas fa-sliders-h","Add Avivator Chart",()=>{
                scShowVivDialog(cm,hconf);

            });
        }
        
        if (cm._ssfilter){
            cm._ssfilter.removeFilter(false);
        }
        let ss= cm.viewData.spatial_stats;
        if (!ss){
            ss={};
            cm.viewData.spatial_stats=ss;
        }
        for (let n in cm.dsIndex){
            const ds =cm.dsIndex[n];
            if (ds.custom.spatial_stats && !(ss[ds.name])){
                ss[ds.name]={show_points:true,show_density:false};
            }
        }
            
    
        
        for (let n in cm.dsIndex){
            const ds= cm.dsIndex[n];
            if (ds.custom.spatial_stats && cm.viewData.initialCharts[ds.name]){
            
                scAddSSMenuIcon(ds,cm)
            
            }
        }

    scUpdateImageConfigs(hconf.images,cm,urlbase);
    }
}

function scAddSSMenuIcon(ds,cm){
    cm.addMenuIcon(ds.name,"fas fa-chart-pie","Add Spatial Chart",()=>{
                    scShowAddSpatialChart(ds,cm);
    });
    cm.addMenuIcon(ds.name,"fas fa-cog","Adjacency Plot Settings",()=>{
                    scShowCellAdjacencyDialog(ds,cm);
    });
    cm.addMenuIcon(ds.name,"fas fa-project-diagram","Show Single Cell cell Interaction",()=>{
                    scShowSingleCellInteraction(ds,cm);
    });
    
    cm.addMenuIcon(ds.name,"fas fa-viruses","Select Cell Types To Display",()=>{
        scShowCellTypesToDisplay(ds,cm);
});
}

function scShowVivDialog(cm,hconf){


const ot = hconf.avivator["ome_tiff_files"];
const url = cm._urlbase+"/" ;
const options=[];
for (let sid in ot){
    options.push({text:sid,value:sid})
}
const conf = {
    title:"Add Channel",
    controls:[
        {
            label:"ROI",
            type:"dropdown",
            id:"dd1",
            items:options
        }
    ],
    buttons:[
        {
            text:"Add",
            method:(c)=>{
                scAddVivChart(cm,hconf,{roi: c["dd1"],image:url+ot[c["dd1"]]});			
            }
        }
        
    ]

}
cm.showCustomDialog(conf);

}

function scAddVivChart(cm,hconf,conf){

const defaults = hconf.avivator.defaults;
const densParam = hconf.default_color_field;

const chart = {
    type:"viv_scatter_plot",
    title:conf.roi,
    background_filter:{
        column:"sample_id",
        category:conf.roi
    },
    color_legend:{
        display:false
    },
    radius:3,
    trim_color_scale:"0.01",
    color_by:densParam,
    param:hconf.position_fields.concat([densParam]),
    roi:{
        max_x:hconf.roi_sizes[conf.roi][0],
        max_y:hconf.roi_sizes[conf.roi][1]
        
    },
    viv:{
        image_properties:defaults,
        url:conf.image
    }
    


}

cm.addChart("cells",chart,false)

}


function scShowSingleCellInteraction(ds,cm){
const cvals =ds.dataStore.getColumnValues("Cell Type 1");

const options = [...cvals].sort().map(x=>{
    return {text:x,value:x};
})
const config ={
    title:"Single Cell Cell Interaction",
    controls:[
        {
            label:"Cell Type 1",
            type:"dropdown",
            id:"ct1",
            items:options
        },
        {
            label:"Cell Type 2",
            type:"dropdown",
            id:"ct2",
            items:options
        }
    ],
    buttons:[
        {
            text:"Change",
            method:(c)=>{
                scShowCellCellInteraction(c["ct1"],c["ct2"],cm,ds);
            
            }
        }
        
    ]
}
cm.showCustomDialog(config);
}

function scShowCellTypesToDisplay(dso,cm){
const ds =dso.dataStore;
const samples = ds.getColumnValues("Cell Type 1");
const celltypes = cm.viewData.spatial_stats[dso.name].celltypes_to_display;
function getCheckBoxes(existing_vals){
    return samples.map(x=>{
        let checked= true;
        if (existing_vals){
            checked = existing_vals.indexOf(x) !== -1;
        }
        return {
            id:x,
            label:x,
            type:"checkbox",
            current_value:checked
        }
    })
}
const cbs1= getCheckBoxes(celltypes?celltypes[0]:null);
const cbs2= getCheckBoxes(celltypes?celltypes[1]:null)

const config= {
    title:"Add Spatial Chart",
    maxHeight:400,
    controls:[
        {
            type:"checkboxgroup",
            label:"Cell Type 1",
            id:"ct1",
            checkboxes:cbs1
        },
        {
            type:"checkboxgroup",
            label:"Cell Type 2",
            id:"ct2",
            checkboxes:cbs2
        }
    ],
    buttons:[
        {
            text:"OK",
            method:(c)=>{
                const arr=[[],[]];
                let index=0;
                for (let n of ["ct1","ct2"]){
                    for (i of samples){
                        if (c[`${n}_${i}`]){
                            arr[index].push(i)
                        }
                    }
                    index++;
                }
                if (arr[0].length === samples.length && arr[1].length === samples.length){
                    delete cm.viewData.spatial_stats[dso.name].celltypes_to_display;
                }
                else{
                    cm.viewData.spatial_stats[dso.name].celltypes_to_display=arr
                }
            }
        }
        
    ]
};
cm.showCustomDialog(config);

}
function scShowAddSpatialChart(dso,cm){
const ds =dso.dataStore;
const s =ds.getColumnList("text");
const cols= s.filter(x=>!x["field"].startsWith("Cell Type")).map(x=>{
    return {text:x.name,value:x.field};
});
const pivot = dso.custom.main_pivot;


const samples = ds.getColumnValues(pivot);
//const charts = ["MoruetaHolme","quadratCounts","PCF_r10","PCF_r20","Network"];
const cs = dso.custom.chart_types.map((x,i)=>{
    return {value:i,text:x.name}
    
});
const dds= samples.map(x=>{
    return {value:x,text:x}
});
const config= {
    title:"Add Spatial Chart",
    controls:[
        {
            label:"Parameter",
            type:"dropdown",
            id:"column",
            items:cols,
            current_value:pivot,
            func:(value,dialog)=>{
                const vs  = ds.getColumnValues(value).map(x=>{
                    return {value:x,text:x}
                });
                dialog.changeDropDownContent("value",vs)
            }

        },
        {
            label:"Show all Values",
            type:"checkbox",
            current_value:false,
            id:"all_values"
        },
        {
            label:"Value",
            type:"dropdown",
            id:"value",
            items:dds
        },
        {
            label:"Chart",
            type:"dropdown",
            id:"chart",
            items:cs
        }
    
    ],
    buttons:[
        {
            text:"Add Chart",
            method:(c)=>{
                const celltypes = cm.viewData.spatial_stats[dso.name].celltypes_to_display;
                const conf = scGetSpatialStatsChart( dso.custom.chart_types[c.chart],c.value,c.column,celltypes,c.all_values)
                cm.addChart(dso.name,conf);
            }
        }
        
    ]
}
cm.showCustomDialog(config);
}


function scGetSpatialStatsChart(chart,sample,column,ct,show_all){
    if (!column){
        column="sample_id";
    }

    const type = chart.type;



    let conf={};
    if (type=="pie"){
            conf =   {
                type:"multi_pie_chart",
                param:["Cell Type 1","Cell Type 2",chart.column,column],
                category:sample,
                title:sample
            }
    }
    else if (type==="heatmap"){
        conf = {
            type:"single_heat_map",
            param:["Cell Type 1", "Cell Type 2",chart.column,column],
            category:sample,
            axis:{
                y:{textSize:13,label:"",size:80},
                x:{textSize:13, label:"",size:80,rotate_labels:true}
            }
            
        };
        if (chart.stat_cutoff){
            conf.param.push(chart.stat_cutoff);
            conf.cutoff=0.05
        }
        if (show_all){
            conf.type="single_heat_map_group";
        }


        /*let cs=[-15,15];
        if (type.startsWith("PCF")){
            cs= [0,2];
        }
        else if (type==="quadratCounts"){
            cs=[-1,1]
        }*/
        conf.color_scale={
                min_max:chart.scale
        }

    }
    else if (type==="network"){
        conf = {
            type:"cell_network_chart",
            stat_cutoff:0.05,
            param:[column,"Cell Type 1", "Cell Type 2",chart.column,chart.stat_cutoff],
            category:sample,   
        };
    
    }
    if (ct){
        conf.specific_only=[ct[0],ct[1]];
    }

    return conf;

}


function scShowCellAdjacencyDialog(ds,cm){
const ss= cm.viewData.spatial_stats[ds.name];
const cds= cm.getDataSource("cells");
//const cat = ds.column_link_to.columns[1].link_to;
let vals= cds.getColumnValues(ds.custom.main_pivot);
vals= vals.map(x=>{
    return {text:x,value:x}
})



const config={
     title:"Adjacency Plot Settings",
     text:"",
    buttons:[{
        text:"Change ROI",
        method:(c)=>{
            scShowSpatialROI(c.roi,cm,ds);
        }
        
    }],
     controls:[
        {
            type:"dropdown",
            label:"ROI to show",
            items:vals,
            id:"roi"
            
        },
        {
             type:"checkbox",
             current_value:ss.show_points,
             label:"Show Points",
             func:(x)=>{
                    const chs =cm.getAllCharts("cells").filter(x=>x.config.id.startsWith("ss_"));
                    ss.show_points=x;
                    for (let ch of chs){
                        const op =x?0.8:0;
                        ch.config.opacity=op
                           ch.app.setPointOpacity(op);
                        ch.app.refresh();
                    }
                }
          },
        {
            type:"checkbox",
            max:10,
            min:0,
            current_value:ss.show_density,
            label:"Show Density Contours",
            func:(x)=>{
                ss.show_density=x;
                const cds = ds.dataStore;
                const hl = cds.getHighlightedData();
                if (! hl){
                    return;
                }
                const info = cds.getRowAsObject(hl[0],["sample_id","Cell Type 1","Cell Type 2"]);
                const chs =cm.getAllCharts("cells").filter(x=>x.config.id.startsWith("ss_"));
                for (let ch of chs){
                    if (!x){
                        ch.removeCategory(1);
                        ch.removeCategory(2);
                    }
                    else{
                        ch.config.category1= info["Cell Type 1"];
                        ch.config.category2= info["Cell Type 2"];
                        ch.onDataFiltered();
                    }
                }
            }
        }
    ]
};

cm.showCustomDialog(config);

}

function scGetSpatialROITables(ds,offset){
const im = ds.custom.pcf_charts;
  const im_chart={
         "type":"image_table_chart",
      "param":[
        im.key_column
      ],
      "images":{
        "base_url":im.base_url,
        "type":"png"
      },
      "id":"ss_image_plot",
      "size":[
        800,
        320
      ],
      "image_width":300,
    "title":"PCF plots",
   
      "position":[
        10,
        offset+20
      ]
}
const table = {
      
      "title":"",
      "legend":"",
      "type":"table_chart",
      "param":ds.custom.default_table_columns,
      "id":"W2LHAW",
      "size":[
        800,
        320
      ],
      "column_widths":{

      },
      "position":[
       820,
       offset+20
      ]
    
}
return [im_chart,table];


}

function scShowCellCellInteraction(c1,c2,cm,ds){
cm.removeAllCharts();
const charts = scGetSpatialROITables(ds,300);
const dgf = ds.custom.default_group_field;
const val_col= ds.custom.default_group_chart_values

if (ds.custom.main_pivot==="sample_id"){
    charts.push({
        type:"single_series_chart",
        categories:[c1,c2],
        param:[dgf,"Cell Type 1","Cell Type 2","sample_id",val_col],
        scale:[0,2],
        position:[10,10],
        title:val_col+" "+c1+" "+c2,
        size:[400,300]
    })
}

const apcf= ds.custom.average_pcf;
let apcf_chart =null;

if (apcf){
    let urlbase = cm._urlbase.replace("images","");
    const index1= ds.dataStore.getColumnValues("Cell Type 1").indexOf(c1);
    const index2= ds.dataStore.getColumnValues("Cell Type 2").indexOf(c2);
    urlbase=urlbase+ds.name+"_pcf_av";

    const ims=[]
    for (item of apcf.conditions){
        ims.push([
            urlbase+"/"+item[0]+"_"+index1+"_"+index2+".png",
            item[1]
        ])

    }
    apcf_chart={
        type:"image_box",
        images_per_row:apcf.images_per_row,
        image_size:apcf.image_size,
        images:ims,
        size:[600,200],
        position:[420,10]
    }
    
}



cm._getColumnsThen(ds.name,["Cell Type 1","Cell Type 2"],()=>{
    ds._ssfilter.filter("filterMultipleCategories",["Cell Type 1","Cell Type 2"],[c1,c2]);

    for (let ch of charts ){
        cm.addChart(ds.name,ch);
    }
    if (apcf){
        cm.addChart(ds.name,apcf_chart)
    }
});

}

function scShowSpatialROI(roi,cm,ds){
const n = ds.dataStore.getColumnValues("Cell Type 1").length;
const hmsize=80+(8*n);
cm.removeAllCharts();
const pivot = ds.custom.main_pivot;
const tables = scGetSpatialROITables(ds,hmsize);
const celltypes = cm.viewData.spatial_stats[ds.name].celltypes_to_display;
//const fe = scGetSpatialStatsChart("Network",roi,"sample_id",celltypes);
for (let i in ds.custom.default_charts){
    const index =  ds.custom.default_charts[i];
    const ch =  scGetSpatialStatsChart(ds.custom.chart_types[index],roi,pivot,celltypes);
    ch.position= [10+(i*(hmsize+10)),10];
    ch.size=[hmsize,hmsize];
    tables.push(ch)
}

/*const mh = scGetSpatialStatsChart("MoruetaHolme",roi,"sample_id",celltypes);
mh.position=[10,10];
mh.size=[hmsize,hmsize];

const qc= scGetSpatialStatsChart("quadratCounts",roi,"sample_id",celltypes);
qc.position=[20+hmsize,10];
qc.size=[hmsize,hmsize];

const pcf= scGetSpatialStatsChart("PCF_r10",roi,"sample_id",celltypes);
pcf.position=[30+hmsize*2,10]
pcf.size=[hmsize,hmsize];
*/

for (let ch of tables){
    cm.addChart(ds.name,ch);
}

cm._getColumnsThen(ds.name,[pivot],()=>{
    ds._ssfilter.filter("filterCategories",[pivot],[roi]);
    if (pivot==="sample_id"){
        const settings = cm.viewData.spatial_stats[ds.name];
        const cid = "ss_"+roi;
        const cat = ds.column_link_to.columns[1].link_to;

        
        const details={
            opacity:settings.show_points?0.8:0,
            cats:null,
            point_size:5
        }
        
        scAddImage(roi,cm._urlbase,cm,cm.hyperion_config,{color_by:cat,id:cid,details:details,size:[500,500]});
        const he = n*8 + 30;
        const rc =   {
            "title":cat,
            "type":"row_chart",
            "param":cat,
            "size":[
                250,
                he
            ],
            "position":[
                520,
                10
            ]
            }
        cm.addChart("cells",rc);
    }
})

}



function scShowCellAbundanceDialog(cm,hyp_conf){
const dn = {};
const ds = cm.getDataSource("cells");
const p1_cols = ds.getColumnList("text").map(x=>{
    return {text:x.name,value:x.field};
});
const p3_cols= ds.getColumnList("text").map(x=>{
    return {text:x.name,value:x.field};
});

const p3_values= ds.getColumnValues(p3_cols[0].value).map(x=>{
    return {text:x,value:x};
})

for (let sample in hyp_conf.roi_sizes){
    const dim = hyp_conf.roi_sizes[sample];
    dn[sample]=(dim[0]*dim[1])/1000000;
}
const conf={
    type:"custom_box_plot",
    denominators:dn,
    axis:{
        y:{label:"cells / mm2"}
    }
}
cm.getDataSource("cells")
const controls= [
    {
        label:"Groups (X axis)",
        type:"dropdown",
        id:"p1",
        items:p1_cols
    },
    {
        label:"Selection Column",
        type:"dropdown",
        id:"p3",
        items:p3_cols,
        func:(value,dialog)=>{
            const vs  = ds.getColumnValues(value).map(x=>{
                return {value:x,text:x}
            });
            dialog.changeDropDownContent("value",vs)
        }

    },
    {
        label:"Value",
        type:"dropdown",
        id:"value",
        items:p3_values
    },
];
cm.showCustomDialog({
    title:"Add Abundance Chart",
    controls:controls,
    buttons:[
        {
            text:"Add Chart",
            method:(c)=>{
                conf.param= [c["p1"],"sample_id",c["p3"]];
                conf.category=c["value"];
                conf.title = c["value"];			
                cm.addChart("cells",conf);
            }
        }
        
    ]
})

}

function scShowBulkChangeImagesDialog(cm,conf){
//get all the image plots
let  plots = cm.getAllCharts("cells");
plots=plots.filter(x=>x.config.background_filter && x.config.type==="density_scatter_plot");
//the defualt value from the first chart
const c= plots[0].config 
//get the columns to color by
let cols = cm.dsIndex["cells"].dataStore.getColumnList();
cols=cols.map(x=>{
    return {
        text:x.name,
        value:x.field
    }
});
cols.push({
    text:"None",
    value:"_none"
});

//get all image names
const is = new Set();
for (let s in conf.images){
    for (let i of  conf.images[s] ){
        is.add(i[0])
    }
}
const ims = [];
for (let i of is){
    ims.push({text:i,value:i})
}


const config={
     title:"Bulk Alter Image Plots",
     text:"",
     controls:[
        {
             type:"slider",
                max:1,
                min:0,
                current_value:c.opacity,
                label:"Point Opacity",
                func:(x)=>{
                    for (let p of plots){
                        p.config.opacity=x;
                        p.app.setPointOpacity(x)
                        p.app.refresh();
                    }
                }
          },
        {
            type:"slider",
            max:10,
            min:0,
            current_value:c.radius,
            label:"Point Size",
            func:(x)=>{
                for (let p of plots){
                    p.config.radius=x;
                    p.app.setPointRadius(x)
                    p.app.refresh();   
                }
            }
        },
        {
            label:"Color By",
            type:"dropdown",
            items:cols,
            current_value:c.color_by || "_none",
            func:(x)=>{
                for (let p of plots){
                      if (x==="_none"){
                            delete p.config.color_by
                            p.colorByDefault();
                        }
                            else{
                            p.config,color_by=x;
                            p.colorByColumn(x);
                        }
                }
            }
        },
        {
            type:"dropdown",
                current_value:c.background_image.name,
                items:ims,
                label:"Change Image",
                func:x=>{
                    for (let p of plots){
                        const cf =p.config;
                        let t  = cf.title.split("-")[0];
                        let f = null;
                        for (let i of cf.image_choices){
                            //get the image of this name
                            if(i[0]===x){
                                f=i;
                                break;
                            }
                        }
                        //image does not exist - ignore
                        if (!f){
                            continue;
                        }
                        cf.background_image.url=f[1];
                        cf.background_image.name=f[0];
                        if (f[2]){
                                cf.background_image.width=f[2];
                                cf.background_image.height=f[3];
                                cf.background_image.position=[f[4],f[5]];
                    }
                    p.setTitle(t+"-"+f[0])
                    p.addBackgroundImage(cf.background_image,true);
                    }
                }
            }
        ]
}
cm.showCustomDialog(config);
}


function scUpdateImageConfigs(images,cm,urlbase){
const charts = cm.getAllCharts("cells");
for (let ch of  charts){
    const c = ch.config;
    if (c.background_filter && c.type==="density_scatter_plot"){
        const ims =images[c.background_filter.category];
        const imc = c.image_choices;
        for (let im of ims){
            if (!imc.find(x=>x[0]===im[0])){
                imc.push([im[0],`${urlbase}/${im[1]}.png`,im[2],im[3],im[4],im[5]]);
            }
        }
    }
}
}



function scShowUploadImageDialog(cm,images,urlbase){
const url = "/meths/execute_project_action/"+project_id;
const data = {
    method:"upload_images",
    args:{}
}
new scFileUploadDialog({
    url:url,
    onupload:()=>data,
    onfinished:response=>{
                const n_i = response.data;
                for (let sid in n_i){
                    let isid= images[sid];
                    if (!isid){
                        images[sid]= n_i[sid]
                    }
                    else{
                        for (let im of n_i[sid]){
                            if (!isid.find(x=> x[0]===im[0])){
                                isid.push(im)
                            }
                        }
                        
                    }
                }
                   scUpdateImageConfigs(images,cm,urlbase)
    }
});
}


function scGetLocalViewLoader(folder){
    return async (name)=>{
      const resp= await fetch(`${folder}/views.json`);
     
      const js=  await resp.json();  
      return js[name];
    }
    
}

function scGetLocalDataLoader(folder,dataSources){
    const loaders = {}
    for (let ds of dataSources){
        loaders[ds.name]=new DataLoader(`${folder}/data/${ds.name}.b`,ds.size)
    }
    return async(columns,dataSource,size)=>{
      return await loaders[dataSource].getColumnData(columns)
     
      
    }
    
}

async function scViewLoader(name){
  const resp= await fetch("/meths/execute_project_action/"+project_id,
    {
        method:"POST",
        body:JSON.stringify({
            method:"get_view",
            args:{
                view:name
            }
            }),
        headers:{
            "Accept":"application/json,text/plain,*/*",
            "Content-Type":"application/json"
        }
    });
    const js=  await resp.json();
    // ***** legacy stuff
    if (!js.data.initialCharts){
        return {
            initialCharts:js.data
        }
    }
    // ********
    return js.data;
}

function scGetDefaultImages(cm){
    const d= {}
    const chs = cm.getAllCharts("cells");
    chs.forEach(x=>{
        const c = x.config
        if (c.background_image){
            d[c.background_filter.category]=c.background_image;
        }
    });
    return d;
}


function scListener(type,cm,data,proj){
    switch(type){
        case "state_saved":
            

            const c = proj.data.hyperion_config;
            data.default_images= c.default_images || {};
            executeProjectAction("save_state",{state:data}).then(resp=>{
                if (resp.success){
                    cm.createInfoAlert("Data Saved",{duration:2000});
                }
                else{
                    cm.createInfoAlert("UnableToSaveData",{duration:3000,type:"danger"});
                }
            })	
            break;
        case "view_loaded":
            scChangeURLParam("view",cm.currentView);
            scInitApp(cm,proj);
            break;
        case "chart_added":
            const chart=data;
            if (chart.config.type.endsWith("_scatter_plot")){
                chart.config.trim_color_scale="0.01";
                const col = proj.data.hyperion_config.default_color_field;
                if (!chart.config.color_by && col){
                    chart.colorByColumn(col);
                }
            }
            if (chart.config.background_image){
                scAddSetDefaultChoice(chart,proj.data.hyperion_config);
            }
    }
}


function scAddSetDefaultChoice(ch,hyp_conf){
 ch.addToContextMenu = ()=>{
  
 return[{
  text:"Set as Default Image",
  icon:"fas fa-file-image",
  func:()=>{
    const c  = ch.config;
   const s = c.background_filter.category;
    hyp_conf.default_images[s]= c.background_image;
  }
}];

}

}


class SCDataUploader{
constructor(){
    this.showUploadDialog();
}
    
uploadData(file,fields,has_headers,delimiter){
    let url = "/meths/execute_project_action/"+project_id;
    let self = this;
    let data = {
            method:"upload_data_file",
            arguments:{
                fields:fields,
                has_headers:has_headers,
                delimiter:delimiter
            }
    }
    this.waiting_icon= new WaitingDialog("Uploading And Processing File");
    this.waiting_icon.wait("Uploading File");
    mlvUploadFile(file,url,data,function(response){		
        self.dataUploaded(response);	
    });
    this.mlv_file_upload.remove();		
}


showUploadDialog(){
    let config = {
             compulsory_fields:{
                     1:{label:"Cell ID",datatype:"text"},
                }
            }
    
    this.mlv_file_upload = new MLVFileUploadDialog(config);
    let self = this;
    this.mlv_file_upload.setUploadCallback(function(file,fields,has_headers,delimiter){
        self.uploadData(file,fields,has_headers,delimiter);
    });
}


dataUploaded(response){
    let data = response.data
    if (response.success){
        location.reload()   	
    }
    else{
        this.waiting_icon.showMessage("There was a problem uploading the data","danger");
    }	
}
}


async function executeProjectAction(action,args){
    if (!args){
        args={}
    }
    let data={
        method:action,
        args:args	
    }
    const resp = await fetch("/meths/execute_project_action/" +project_id,
    {
        method: "POST",
        body: JSON.stringify(data),
        headers: {
            "Accept": "application/json,text/plain,*/*",
            "Content-Type": "application/json"
        }
    });
return await resp.json();
}

function scChangeURLParam(param,val){

const url = new URL(window.location);
(url.searchParams.has(param) ? url.searchParams.set(param, val) : url.searchParams.append(param, val));

url.search = url.searchParams;
url == url.toString();
history.pushState({}, null, url);

}


class DataLoader{
    constructor(url,size){
      this.url=url;
      this.index=null;
      this.size=size;
    }
    async getColumnData(cols){
      if (!this.index){
        const iurl = this.url.replace(".b",'.json')
        const resp = await fetch(iurl);
        this.index = await resp.json();
      }
      
      return await Promise.all(cols.map(async (c)=>{
        const i = this.index[c.field];
  
        const resp =await fetch(this.url,
            { headers:
              {
                responseType: "arraybuffer",
                range:`bytes=${i[0]}-${i[1]}`
              }
            });
        const bytes = await resp.arrayBuffer();
        const output = pako.inflate(bytes);
        const sb = new SharedArrayBuffer(output.length)
        const f  = new Uint8Array(sb);
        f.set(output,0);
        return {data:sb,field:c.field};
      }));
    }
}
export {scLoadData,scLoadMetaData};
