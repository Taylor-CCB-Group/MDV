import { createEl, makeDraggable, makeResizable,removeDraggable,removeResizable,createMenuIcon,splitPane} from "../utilities/Elements";
import { BaseChart} from "./BaseChart.js";
import { PopOutWindow } from "../utilities/PopOutWindow";
import  DataStore from "../datastore/DataStore.js";
import  "./HistogramChart.js";
import  "./RowChart.js";
import "./TableChart.js";
import "./WGL3DScatterPlot.js";
import "./WGLScatterPlot.js";
import "./RingChart.js";
import "./TextBoxChart.js";

import {BaseDialog} from "../utilities/Dialog.js";
import {getRandomString} from "../utilities/Utilities.js";
import {csv,tsv,json} from"d3-fetch";


/**
* The object to manage charts
* @param {string|DOMelement} div - The DOM element or id of the element to house the app
* @param {object[]} datasources - An array of datasources, each should have the following
* <ul>
* <li> name  -the name of the datasource</li>
* <li> size - the size (number or rows) of the data set </li>
* <li> columns - a list of objects describing each column see {@link DataStore#addColumn} <li>
* </ul>
* @param {function | object} - either a function describing how to load columns or an object(s) 
* describing where to get the data from
* <ul>
* <li> url - the url of the file <li>
* <li> type - either csv,tsv or json <li>
* <li> dataSource - the name of th datasource <li>
* </ul>
* 
*/
class ChartManager{

    //
    constructor(div,dataSources,dataLoader,config={}){

        this.infoAlerts={};
        

        //each entry in dataSources will contain
        //  dataSource - the actual dataStore object
        //  name - the name given to this data source
        //  menuBar the dom menu associated with this element
        //  contentDiv the div that the charts associated with the datastore willbe added
        this.dataSources=[];
        this.dsIndex={};
        for (const d of dataSources){
            const ds= {
                name:d.name,
                dataStore:new DataStore(d.size,{columns:d.columns})
            }
            this.dataSources.push(ds);
            this.dsIndex[d.name]=ds
        }

        this.transactions={};
       
     
        this.containerDiv= typeof div === "string"?document.getElementById(div):div;
        this.listeners={};

       

        
     

        //each entry in charts will contain
        //  chart - the actual chart
        //  win - the popout wimdow it is in (or null)
        //  dataSource - the data source associated with it 
        this.charts={};


        this.config=config;
        const c = this.config
        c.chart_color=c.chart_color || "white";

        //create menu bar and containers for each datasource
    
        const panes = splitPane(this.containerDiv,{number:dataSources.length});
        for (let n=0;n<this.dataSources.length;n++){
            const p = panes[n];
            p.style.display="flex";
            p.style.flexDirection="column";
            const ds= this.dataSources[n];
            ds.charts=[];
            ds.menuBar = createEl("div",{
                classes:["ciview-menu-bar"]          
            },p);
            this._setUpMenu(ds);
            ds.contentDiv=createEl("div",{
                styles:{
                    flex:"1 1 auto",
                    position:"relative",
                    overflow:"auto",
                    background:dataSources[n].color || "beige"
                }
            },p);
         

        }
        config.initialCharts= config.initialCharts || {}
        if (typeof dataLoader === "function" ){
            this.dataLoader = dataLoader;
            const ic= config.initialColumns;
            if (ic){ 
                for (let dataSource in ic){     
                    this.loadInitialColumns(ic[dataSource],dataSource, config.initialCharts[dataSource] || [])
                } 
            }
        }
        else{
            if (!Array.isArray(dataLoader)){
                dataLoader=[dataLoader]
            }
            for (let item of dataLoader){
                this.loadFile(item,config.initialCharts)
            }
        }
  
    }

    loadFile(info,charts){
        const meths = {csv:csv,json:json,tsv:tsv}
        const iid  =  this.createInfoAlert("loading file",{spinner:true})
        meths[info.type](info.url).then(data=>{
            const cols={};
            const dataSource= info.dataSource
            const ds =this.dsIndex[dataSource].dataStore;
            const all_cols =  ds.getAllColumns();
            for (let f of all_cols){
                cols[f]=[];
            }
            for (let i=0;i<data.length;i++){
                const row =data[i];
                if (i+1%100===0){
                    this.updateInfoAlert(iid,`processed ${i}/${data.length} rows`)
                }
                for (let col in cols){
                    cols[col].push(row[col])
                }
            }
            let proc=0;
            for  (let col in cols){
                ds.setColumnData(col,cols[col]);
                proc++;
                this.updateInfoAlert(iid,`processed ${proc}/${all_cols.length} columns`)

               
            }
            this.updateInfoAlert(iid,"complete",{duration:2000})
            if (charts[dataSource]){
                for (let ch of charts[dataSource]){
                    this.addChart(dataSource,ch)
                }
            }
        })

    }




    loadInitialColumns(columns,dataSource,charts){
        this.loadColumnSet(columns,dataSource,()=>{     
                for (let chart of charts){
                    this.addChart(dataSource,chart);
                }
            
        });
    }


    getState(){
        const initialColumns={};
        const initialCharts={};
        for (const ds of this.dataSources){  
            initialColumns[ds.name]= ds.dataStore.getLoadedColumns();
            initialCharts[ds.name]=[];
        }
        for (let chid in this.charts){
            const chInfo = this.charts[chid];
           
            const chart = chInfo.chart;
            const config = chart.getConfig();
            const div =  chart.getDiv();
            config.position = [div.offsetLeft,div.offsetTop];
        
            initialCharts[chInfo.dataSource.name].push(config);
            
        }
        return{
            initialColumns:initialColumns,
            initialCharts:initialCharts
        }
    }

    createInfoAlert(msg,config={}){
        let id = getRandomString();
        const len = Object.keys(this.infoAlerts).length;
        config.type= config.type || "info"
        const div = createEl("div",{
            classes:["ciview-info-alert","ciview-alert-"+config.type],
            styles:{
                right:"10px",
                top:50+(len*40)+"px",
            },
          
        },this.containerDiv);
        let spinner = null;
        const text=  createEl("span",{text:msg},div);
        if (config.spinner){       
            spinner=createEl("i",{
                classes:["fas","fa-spinner","fa-spin","ciview-info-alert-spin"]
            },div);
        }     
        this.infoAlerts[id]={
            div:div,
            text:text,
            spinner:spinner,
            type:config.type 
        };
        if (config.duration){
            this.removeInfoAlert(id,config.duration)
        }
        return id;
    }

    updateInfoAlert(id,msg,config={}){
        const al =this.infoAlerts[id];
        if (config.type && al.type !==config.type){
            al.div.classList.remove("ciview-alert-"+al.type);
            al.div.classList.add("ciview-alert-"+config.type);
            al.type=config.type;
        }
        al.text.textContent=msg;
        if (config.duration){
            this.removeInfoAlert(id,config.duration);
        }
    }

    removeInfoAlert(id,delay=2000){
        const spinner = this.infoAlerts[id].spinner;
        if (spinner){
            spinner.remove();
        }
        setTimeout(()=>{
            this.infoAlerts[id].div.remove();
            delete this.infoAlerts[id];
            let top =50;
            for (let i in this.infoAlerts){
                this.infoAlerts[i].div.style.top = top+"px";
                top+=40;
            }
        },delay);
    }


    //will load the columns and when loaded will execute the callback,
    //passing a list of any columns which failed to load
    loadColumnSet(columns,dataSource,callback,split=10,threads=2){
        const id = getRandomString();
      
     
        this.transactions[id]={
            callback:callback,
            columns:[],
            totalColumns:columns.length,
            failedColumns:[],
            nextColumn:0,
            columnsLoaded:0,
            id:id
        }
        let col_list=[];
        const t  = this.transactions[id]; 
        for (let col of columns){
            col_list.push(col);
            if (col_list.length===split){
                t.columns.push(col_list);
                col_list=[];
            }
        }
        if (col_list.length!==0){
            t.columns.push(col_list);
            col_list=[];
        }

        t.alertID= this.createInfoAlert(`Loading Columns:0/${columns.length}`,{spinner:true});
        
           
      
        const max = Math.min(t.columns.length,threads);
       
        for (let n=0;n<max;n++){
            this._loadColumnData(t,dataSource)
        }

        
    }

    loadChartSet(chartSet){
        for (const ds in chartSet){
            for (const ch of chartSet[ds]){
                this.addChart(ds,ch)
            }
        }
    }

    _loadColumnData(trans,dataSource){
        const dataStore=  this.dsIndex[dataSource].dataStore;
       
        const col_list = trans.columns[trans.nextColumn++];
        const columns=[];
        for (let col of col_list){
           columns.push(dataStore.getColumnInfo(col));
        }
        columns.sort((a,b)=>{
            const c= (a.datatype==="double" || a.datatype==="integer")?0:1;
            const d= (b.datatype==="double" || b.datatype==="integer")?0:1;
            return c-d;

        })
       
   
        
        
        this.dataLoader(columns,dataSource,dataStore.size).then(resp=>{
            for (let col of resp){
                dataStore.setColumnData(col.field,col.data);
            }
            trans.columnsLoaded++;
        }).catch(error=>{
            console.log(error);
            trans.columnsLoaded++;
            trans.failedColumns.push(columns);
          
        }).finally(()=>{
            const total = trans.columns.length;
            const loaded = trans.columnsLoaded;
            let all_loaded= loaded*col_list.length;
            all_loaded = all_loaded>trans.totalColumns?trans.totalColumns:all_loaded;
            this.updateInfoAlert(trans.alertID,`Loading Columns:${all_loaded}/${trans.totalColumns}`);
            if (loaded>=total){
                this.updateInfoAlert(trans.alertID,"All columns loaded",{duration:2000});
                trans.callback(trans.failedColumns);     
                delete this.transactions[trans.id];
            }
            if (trans.nextColumn<total){
                this._loadColumnData(trans,dataSource)
            }          
        })
    }


    _setUpMenu(ds){
        createMenuIcon("fas fa-chart-bar",{
            tooltip:{
                text:"Add Chart",
                position:"bottom-right"
            },
            func:()=>{
                new AddChartDialog({},{
                    dataStore:ds.dataStore,
                    callback:(config)=>{
                        this.addChart(ds.name,config,true);
                    }
                });
            }
            },ds.menuBar
        );

        createMenuIcon("fas fa-sync-alt",{
            tooltip:{
                text:"Reset All Filters",
                position:"bottom-right"
            },
            func:()=>{
               ds.dataStore.removeAllFilters();
            }
            },ds.menuBar
        );

        if (this.config.permission==="edit"){
            createMenuIcon("fas fa-save",{
                tooltip:{
                    text:"Save",
                    position:"bottom-right"
                },
                func:()=>{
                    const state = this.getState();
                    this._callListeners("state_saved",state)
                }
    
            },ds.menuBar);
        }
      

    }
    addListener(id,func){
        this.listeners[id]=func;
    }
    removeListener(id){
        delete this.listeners[id];
    }

    _callListeners(type,data){
        for (let id in this.listeners){
            this.listeners[id](type,data);
        }
    }

 

    addChart(dataSource,config,notify=false){
        //**convert legacy data*********** 
        if (config.location){
            const l = config.location;
            const b=5;
            config.size=[l.width*90 + l.width*b -b,l.height*40 + l.height*b -b];
            config.position=[(l.x+1)*b + l.x*90, (l.y+1)*b + l.y*40];
        }
         //**convert legacy data***********

        const ds  = this.dsIndex[dataSource];
        let width=300,height= 300;
        let left=10,top=10;
        if (config.size){
            width=config.size[0];
            height=config.size[1];
        }
        if (config.position){
            left=config.position[0];
            top=config.position[1];
        }

        const chartType= BaseChart.types[config.type];
        if (!chartType){
           // throw `${config.type} chart type does not exist`
          console.log(`${config.type} chart type does not exist`);
          return;
        }
        const div= createEl("div",{
            styles:{
                position:"absolute",
                width:width+"px",
                height:height+"px",
                left:left+"px",
                top:top+"px",
                background:this.config.chart_color
            }
        },ds.contentDiv);
        
       
        
        const chart = new chartType.class(ds.dataStore,div,config);
        this.charts[chart.config.id]={
            chart:chart,
            dataSource:ds
        }


        this._makeChartRD(chart);
        chart.popoutIcon = chart.addMenuIcon("fas fa-external-link-alt","popout",{
            func:()=>{
                this._popOutChart(chart,ds.contentDiv);
            }
        });

        
        chart.addMenuIcon("fas fa-trash","remove chart")
            .addEventListener("click",()=>{   
                chart.remove();
                div.remove();
                delete this.charts[chart.id];
                this._callListeners("chart_removed",chart);
            });


        if (notify){
            this._callListeners("chart_added",chart);
        }
        
        return chart;
      
    }

    getChart(id){
        const cinfo = this.charts[id];
        if (!cinfo){
            return null;
        }
        return cinfo.chart;
    }

    setChartsAsGrid(rowLength=5,size=[300,300],margin=10){
        let top=margin;
        let left =margin;
        let rowSize=0;
        for (let id in this.charts){
            const info = this.charts[id];
            const d= info.chart.getDiv();
            d.style.left=left+"px";
            d.style.top=top+"px";
            //info.chart.setSize(size[0],size[1]);
            left+=size[0]+margin;
            rowSize++;
            if (rowSize===rowLength){
                rowSize=0;
                left=margin;
                top+=size[1]+margin;

            }

            
        }
    }

    addButton(text,callback,tooltip){
        createEl("span",{
            classes:["ciview-button"],
            text:text,
            styles:{
                position:"fixed",
                bottom:"40px",
                right:"40px",
                fontSize:"18px",
                zIndex:100
            }
        },this.containerDiv)
        .addEventListener("click",()=>callback())
    }

    _popOutChart(chart){
        const div= chart.getDiv();
        const chInfo= this.charts[chart.config.id];
        const details={dim:[chart.config.size[0],chart.config.size[1]],pos:[div.style.left,div.style.top]};
        removeResizable(div);
        removeDraggable(div);
        const win = new PopOutWindow(
            //new window opens
            (doc,box)=>{
           
              chart.setSize(box.width,box.height);
              div.style.top="5px";
              div.style.left="5px";
              doc.body.append(div)
              chart.changeBaseDocument(doc)
              doc.body.style.overflow="hidden";
              chart.popoutIcon.style.display="none";
        
            },
            //new window closes
            (doc,box)=>{
              chInfo.dataSource.contentDiv.append(div)
              chart.changeBaseDocument(document);
              div.style.left = details.pos[0];
              div.style.top= details.pos[1];
              chart.setSize(details.dim[0],details.dim[1]);
              this._makeChartRD(chart);
              chart.popoutIcon.style.display="inline";
              delete chInfo.window
              
            },
            //config
            { 
                onresize:(doc,box)=>{
                    chart.setSize(box.width,box.height)
                },
                url:"/"
        
            }
        );
        chInfo.window=win;

    }

    _makeChartRD(chart){
        const div = chart.getDiv();
        makeDraggable(div,{
            handle:".ciview-chart-title",
            contain:"topleft"
            //contain:true
        });
        makeResizable(div,{
            onresizeend:(width,height)=>chart.setSize(width,height)
        })
      
    }
}

class ChooseColumnDialog extends BaseDialog{
    constructor(config,content){
        config.footer=true;
        config.width=380;
        super(config,content);
    }
    init(content){}
}

class AddChartDialog extends BaseDialog{
    constructor(config,content){
        config.title= "Add Chart";
        config.columns=2;
        config.footer=true;
        config.width=380;
        super(config,content);
    }
    init(content){
        const types=[];
        for (let type in BaseChart.types){
            const t = BaseChart.types[type]
            types.push({
                name:t.name,
                type:type,

            });
        }
        this.dataStore= content.dataStore;
        types.sort((a,b)=>a.name.localeCompare(b.name));
        this.defaultType=types[0].type;

        createEl("div",{
            text:"Chart Type",
            classes:["ciview-title-div"]
        },this.columns[0]);

        this.chartType = createEl("select",{
            styles:{
                maxWidth:"200px"
            }
        });
        for (let item of types){
            createEl("option",{
                text:item.name,
                value:item.type
            },this.chartType)
        }
        createEl("div",{},this.columns[0]).append(this.chartType);
        this.chartType.addEventListener("change",(e)=>{
            this.setParamDiv(this.chartType.value,content.dataStore);
        });

        createEl("div",{
            text:"Title",
            classes:["ciview-title-div"]
        },this.columns[0]);

        this.chartName= createEl("input",{styles:{width:"150px"}},this.columns[0]);

        createEl("div",{
            text:"Description",
            classes:["ciview-title-div"]
        },this.columns[0]);
        this.chartDescription= createEl("textarea",{styles:{width:"150px",height:"100px"}},this.columns[0]);
      

        createEl("div",{
            text:"Columns",
            classes:["ciview-title-div"]
        },this.columns[1]);
        this.paramDiv = createEl("div",{},this.columns[1]);
        this.setParamDiv(types[0].type,content.dataStore);


        createEl("span",{
            text:"Add",
            classes:["ciview-button"]
        },this.footer).addEventListener("click",()=>this.submit(content.callback));

    }

    submit(callback){
        const config={
            title:this.chartName.value,
            legend:this.chartDescription.value,
            type:this.chartType.value,
            param:this.paramSelects.map((x)=>x.value)
        }
        callback(config);
        this.chartName.value="";
        this.chartDescription.value="";
        this.chartType.value= this.defaultType;
        this.setParamDiv(this.defaultType)

    }

    setParamDiv(type){
        this.paramDiv.innerHTML="";
        const params = BaseChart.types[type].params;
        this.paramSelects=[];
        for (let p of params){
            const d = createEl("div",{styles:{padding:"4px"}},this.paramDiv)
            const sp =createEl("div",{text:p.name+":"},d);
            const holder =createEl("div",{},this.paramDiv);
            const dd = createEl("select",{
                styles:{
                    maxWidth:"200px"
                }
            },holder);
            const ps= this.dataStore.getColumnList(p.type);
            for (let item of ps){
                createEl("option",{text:item.name,value:item.field},dd)
            }
            this.paramSelects.push(dd)
             
        } 
    }
}





export default ChartManager;

export {AddChartDialog};