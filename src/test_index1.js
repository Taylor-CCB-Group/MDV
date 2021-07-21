import "microtip/microtip.min.css";
import "./css/fontawesome-5.15.3/all.min.css";
import 'nouislider/dist/nouislider.min.css'
import "./utilities/css/ContextMenu.css";
import "./charts/css/charts.css";
import "./webgl/css/wgl2di.css";
import "./table/css/slickgrid.css";

import ChartManager from './charts/ChartManager.js';

//
const dataSources1=[
  {
    name:"Cells",
    size:209979,
    columns:[
      {
        name:"PC1",
        field:"pc1",
        datatype:"double"
      },
      {
        name:"PC2",
        field:"pc2",
        datatype:"double"
      },
      {
        name:"PC3",
        field:"pc3",
        datatype:"double"
      },
      {
        name:"PC4",
        field:"pc4",
        datatype:"double"
      },
      {
        name:"PC5",
        field:"pc5",
        datatype:"double"
      },
      {
        name:"PC6",
        field:"pc6",
        datatype:"double"
      },
      {
        name:"PC7",
        field:"pc7",
        datatype:"double"
      },
      {
        name:"PC8",
        field:"pc8",
        datatype:"double"
      },
      {
        name:"names",
        field:"names",
        datatype:"unique",
        stringLength:30
      }
    ]
  },
  {
    name:"Locations",
    size:209979,
    columns:[
      {
        name:"PC1",
        field:"pc1",
        datatype:"double"
      },
      {
        name:"PC2",
        field:"pc2",
        datatype:"double"
      },
      {
        name:"PC3",
        field:"pc3",
        datatype:"double"
      },
      {
        name:"PC4",
        field:"pc4",
        datatype:"double"
      },
      {
        name:"PC5",
        field:"pc5",
        datatype:"double"
      },
      {
        name:"PC6",
        field:"pc6",
        datatype:"double"
      },
      {
        name:"PC7",
        field:"pc7",
        datatype:"double"
      },
      {
        name:"PC8",
        field:"pc8",
        datatype:"double"
      }
      

    ]
  }
];

const config1={
  initialColumns:{
    "Cells":["pc1","pc2","names"],
  },
  initialCharts:{
    "Cells":[
      {type:"bar_chart",param:"pc1"}
    ]
  },
  permission:"edit"
}

function dataLoader1(columns,dataSource,size){

    return new Promise((resolve,reject)=>{
      fetch(`tests/strings.b`)
      .then(response=>response.arrayBuffer())
      .then(data=>{
        const dataList= [];
        let offset=0;
        let arr_len=size;
        for (let column of columns){
          let arrayType=Float32Array;
          let bytes=4;
          if (column.datatype=="unique" || column.datatype=="text"){
              arrayType=Uint8Array;
              bytes=1;
              if(column.datatype=="unique"){
                bytes=column.stringLength;
                arr_len=size*column.stringLength;
              }
             
          }
          const len  = size*bytes;     
          const arr = new arrayType(data,offset,arr_len);     
          const sab = new SharedArrayBuffer(len);
          const new_arr =  new arrayType(sab)
          const t= performance.now();
          for (let n=0;n<len;n++){
              new_arr[n]=arr[n];
          }
          dataList.push({data:sab,field:column.field})
          offset+=len;
        }       
        resolve(dataList);
      }).catch(reject);   
    });

}

const dataSources=[{
  size:384,
  name:"cells",
  columns:[
  {
      "name":"cells",
      "datatype":"unique"
   },
     {
      "name":"Sample",
      "datatype":"text"
   },
        {
      "name":"cell_types",
      "datatype":"text"
   },
           {
      "name":"highLow_runs",
      "datatype":"text"
   },
  {
      "name":"IterativeLSI_UMAP1",
      "datatype":"double"
   },
  {
      "name":"IterativeLSI_UMAP2",
      "datatype":"double"
   },
   {
      "name":"Clusters",
      "datatype":"text"
   }]
}];


const dataLoader = {
  url:"/_tests/testing.umap.cellData.csv",
  dataSource:"cells",
  type:"csv"
} 

const config={
  initialCharts:{
    "cells":[
      {
        type:"wgl_scatter_plot",
        param:["IterativeLSI_UMAP1","IterativeLSI_UMAP2"],
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
      },
      { 
        
        type:"text_box_chart",
        title:"What to Do",
        size:[200,200],
        position:[420,10],
        text:"select clusters and.....",size:[200,200]
      }
    ]
  }
}

const cm = new ChartManager("app1",dataSources,dataLoader,config);
let i =1;
cm.addButton("Proceeed",()=>{
  const ch = cm.getChart("my-clusters");
  const clusters = ch.getFilter();
  if (clusters.length==0 || clusters.length>3){
    cm.createInfoAlert("Please select between 1 and 3 clusters",{type:"warning",duration:2000})
  }
  else{
    //do something with the clusters
  }
 
});







cm.addListener("udu",(type,data)=>{
  switch(type){
    case "state_saved":
      console.log(data);
      cm.createInfoAlert("All saved:"+i++,{duration:2000,type:"danger"})
      break;

  }

});


