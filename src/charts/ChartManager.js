import { createEl, makeDraggable, makeResizable,MDVProgress,removeDraggable,removeResizable,createMenuIcon,splitPane} from "../utilities/Elements";
import  BaseChart from "./BaseChart.js";
import { PopOutWindow } from "../utilities/PopOutWindow";
import  DataStore from "../datastore/DataStore.js";
import CustomDialog from "./dialogs/CustomDialog.js";
import { ContextMenu } from "../utilities/ContextMenu";
import  "./HistogramChart.js";
import  "./RowChart.js";
import "./TableChart.js";
import "./WGL3DScatterPlot.js";
import "./WGLScatterPlot.js";
import "./RingChart.js";
import "./TextBoxChart.js";
import "./HeatMap.js";
import "./ViolinPlot.js";
import "./BoxPlot.js";
import "./SankeyChart.js";
import "./MultiLineChart.js";
import "./DensityScatterPlot";



import {BaseDialog} from "../utilities/Dialog.js";
import {getRandomString} from "../utilities/Utilities.js";
import {csv,tsv,json} from"d3-fetch";
import LinkDataDialog from "./dialogs/LinkDataDialog.js";
import AddColumnsFromRowsDialog from "./dialogs/AddColumnsFromRowsDialog.js";
import ColorChooser from "./dialogs/ColorChooser";



const themes={
    "Dark":{
        title_bar_color:"#222",
        main_panel_color:"black",
        text_color:"white",
        background_color:"#333"
    },
    "Light":{
        title_bar_color:"white",
        main_panel_color:"#f1f1f1",
        text_color:"black",
        background_color:"#bababa"

    }

}



/**
* The object to manage charts
* @param {string|DOMelement} div - The DOM element or id of the element to house the app
* @param {Object[]} datasources - An array of datasources, each should have the following
* <ul>
* <li> name  -the name of the datasource</li>
* <li> size - the size (number or rows) of the data set </li>
* <li> columns - a list of objects describing each column see {@link DataStore#addColumn} <li>
* <li> columnGroups - A list of objects with 'name' and 'columns' field
* </ul>
* @param {function | Object} dataloader - either a function describing how to load columns or an object
* (list of objects) describing the files containing the daata
* <ul>
* <li> url - the url of the file <li>
* <li> type - either csv,tsv or json <li>
* <li> dataSource - the name of the datasource <li>
* </ul>
* A dataLoader can be included as the first item in the list, if dynamic loading as 
* well as static loading 
*
* @param {Object} config extra settings
* @param {Object[]} [config.initialCharts] A list of chart configs to
* initially load
* @param {function} [config.metaDataLoader] A method which returns a Promise
* with a list of the column metadata, given the name of the datasource
* and a list of column field/ids
* @param {function} [config.onDataLoaded] A method which is called once all
* the data is loaded (passing this as the only parameter)
* 
*/
class ChartManager{

    //
    constructor(div,dataSources,dataLoader,config={},listener=null){
        this.listeners={};
        this.infoAlerts={};
        this.progressBars={};
        this.theme="Light"
        // each entry in dataSources will contain
        //  dataSource - the actual dataStore object
        //  name - the name given to this data source
        //  menuBar the dom menu associated with this element
        //  contentDiv the div that the charts associated with the datastore will be added
        this.dataSources=[];
        this.dsIndex={};
        this.columnsLoading={};
        for (const d of dataSources){
            const ds= {
                name:d.name,
                dataStore:new DataStore(d.size,{columns:d.columns,columnGroups:d.columnGroups,offsets:d.offsets,links:d.links}),
                link_to:d.link_to,
                index_link_to:d.index_link_to,
                color:d.color || themes[this.theme].background_color,
                images:d.images,
                genome_browser:d.genome_browser,
                column_link_to:d.column_link_to,
                links:d.links,
                custom:d.custom || {}
            }
            this.dataSources.push(ds);
            this.dsIndex[d.name]=ds;
            this._addDSListeners(ds);
            this.columnsLoading[d.name]={};
            
        }

        if (listener){
            this.addListener("_default",listener)
        }

        this.transactions={};
       
        
        //set up container and top(main menu)
        this.containerDiv= typeof div === "string"?document.getElementById(div):div;
        this.containerDiv.style.display="flex";
        this.containerDiv.style.flexDirection="column";
        
        this.menuBar = createEl("div",{
            classes:["ciview-main-menu-bar"]          
        },this.containerDiv);

        this.leftMenuBar= createEl("span",{},this.menuBar);

        if (config.permission==="edit"){
            createMenuIcon("fas fa-save",{
                tooltip:{
                    text:"Save",
                    position:"bottom-right"
                },
                func:()=>{
                    const state = this.getState();
                    this._callListeners("state_saved",state)
                }
    
            },this.leftMenuBar);
        }

        if (config.all_views){
            this.viewSelect = createEl("select",{},this.menuBar);
            for (let v of config.all_views){
                createEl("option",{text:v,value:v},this.viewSelect)
            }
            this.viewSelect.addEventListener("change",(e)=>{
                this.showSaveViewDialog(()=>this.changeView(this.viewSelect.value));
            })
        }

        if (config.permission==="edit" && config.all_views){
            createMenuIcon("fas fa-plus",{
                tooltip:{
                    text:"Create New View",
                    position:"bottom-right"
                },
                func:()=>{
                    this.showSaveViewDialog(()=>this.showAddViewDialog());
                }
    
            },this.menuBar);
        }

        createMenuIcon("fas fa-palette",{
            tooltip:{
                text:"Change Color Scheme",
                position:"bottom-right"
            },
            func:()=>{
                new ColorChooser(this);
            }

        },this.menuBar);

        
        createMenuIcon("fas fa-adjust",{
            tooltip:{
                text:"Change Theme",
                position:"bottom-right"
            },
            func:(e)=>{
                this.themeMenu.show(e);
            }

        },this.menuBar);

        this._setupThemeContextMenu();
      
        this.contentDiv=createEl("div",{
            styles:{
                flex:"1 1 auto",
                position:"relative"
            }
        },this.containerDiv);
 

     

        //each entry in charts will contain
        //  chart - the actual chart
        //  win - the popout wimdow it is in (or null)
        //  dataSource - the data source associated with it 
        this.charts={};


        this.config=config;
        const c = this.config
        c.chart_color=c.chart_color || "white";

        //load any files first

        this.dataLoader = dataLoader.function;
        this.viewLoader = dataLoader.viewLoader;

        if (dataLoader.files){     
            this.filesToLoad=dataLoader.files.length;
            for (let item of dataLoader.files){
                this.loadFile(item,()=>{
                   this.filesToLoad--;
                   if (this.filesToLoad===0){
                    this._loadView(config,dataLoader,true)
                   }
                });
            }
        }
        else{
            this._loadView(config,dataLoader,true);
        }
         
    }

    _setupThemeContextMenu(){

        this.themeMenu = new ContextMenu(()=>{
            const mItems=[];
           for (let t in themes){
                mItems.push(this.__getMenuItem(t))

           }
           return mItems;
        
        })

    }
    __getMenuItem(theme){
        return {
            text:theme,
            ghosted:this.theme===theme,
            func:()=>this.setTheme(theme)
        }
    }

    setTheme(theme){
        this.theme=theme;
        for (let d of this.dataSources){
            if (d.contentDiv){
                d.contentDiv.style.background= themes[theme].background_color;
            }
            
        }
        for (let id in this.charts){
            this._setChartTheme(this.charts[id].chart)
        }
    }


   

    


    _setChartTheme(chart){
        const t = themes[this.theme];
      
        chart.contentDiv.style.background=t.main_panel_color;
        chart.contentDiv.style.color=t.text_color;
        chart.titleBar.style.background=t.title_bar_color;
        chart.titleBar.style.color=t.text_color;
        if (chart.themeChanged){
            chart.themeChanged();
        }
    }

    _sync_colors(from,to){
        const columns = to.column_link_to.columns;
        for (let item of columns){
            const from_col=from.dataStore.columnIndex[item.link_to];
            const to_col = to.dataStore.columnIndex[item.col];
            const newColors = new Array(to_col.values)
            const colors = from.dataStore.getColumnColors(item.link_to)
            for (let i=0;i<from_col.values.length;i++){
                const val = from_col.values[i];
                const index = to_col.values.indexOf(val);
                if (index!==-1){
                    newColors[index]=colors[i]
                }
            }
            to_col.colors=newColors;

        }
    }

    _initiateOffsets(dataSource){
        const ds = dataSource.dataStore
        const o = ds.offsets;
        const p = o.param;
        const cols= [p[0],p[1],o.groups];
        if (o.background_filter){
            cols.push(o.background_filter);
        }
        this._getColumnsThen(dataSource.name,cols,()=>{
            ds.initiateOffsets();
            ds.updateColumnOffsets();
            ds.updateColumnOffsets(null,true,true);
        })
       

    }

    //load the view metadata or use initialCharts then call _init to load the view 
    _loadView(config,dataLoader,firstTime=false){
       
        //load view, then initialize
        if (config.initial_view){
            this.currentView=config.initial_view;
            this.viewSelect.value = config.initial_view;
            dataLoader.viewLoader(config.initial_view).then(data=>{
                this._init(data,firstTime);
            })     
        }
        //only one view hard coded in config
        else{
            this._init(config.only_view,firstTime)    
        }
    }

    getDataSource(name){
        return this.dsIndex[name].dataStore;
    }

    _init(view,firstTime=false){

        const dsToView=[];
        let charts=[];
        //no initial view
        if (!view){
            for (let ds in this.dsIndex){
                dsToView.push(ds);
            }
            this.viewData={};
        }

        else{
            this.viewData= view;
            charts= view.initialCharts
            for (let ds in charts){
                dsToView.push(ds);
            }
        }

        for (let ds of this.dataSources){
            delete ds.contentDiv;
            delete ds.menuBar;
            if (ds.column_link_to){
                this._sync_colors(this.dsIndex[ds.column_link_to.dataSource],ds);
            }
        }
        const col = themes[this.theme].background_color;
        //add all the appropriate panes (one per datasource)
        const panes = splitPane(this.contentDiv,{number:dsToView.length});
        for (let n=0;n<dsToView.length;n++){        
            const p = panes[n];
            p.style.display="flex";
            p.style.flexDirection="column";
            const ds= this.dsIndex[dsToView[n]];
            this.columnsLoading[ds.name]={};
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
                    background:col
                }
            },p);
        }

        
         //any first time initiation
         if (firstTime){
            //initiate offsets if any
            for (let d of this.dataSources){
                if (d.dataStore.offsets){
                   this._initiateOffsets(d)
                }
            }

        }


        //need to create a set to create track of 
        //charts loaded
        this._toLoadCharts = new Set();
        for (let ds in charts){         
            for (let ch of charts[ds]){
                this._toLoadCharts.add(ch);
            }
        }
        //nothing to load - call any listeners
        if (this._toLoadCharts.size==0){
            delete this._toLoadCharts;
            this._callListeners("view_loaded",this.currentView)
        }
        //add charts - any columns will be added dynamically
        for (let ds in charts){  
            for (let ch of charts[ds]){
                this.addChart(ds,ch);                            
            }
        }

    }

    _addDSListeners(ds){
        ds.dataStore.addListener("l1",(type,data)=>{
            if (type==="column_removed"){
                this._columnRemoved(ds,data)
            }
            else if (type ==="data_highlighted"){
                this._callListeners(type,data);
            }
            else if (type==="filtered"){
                if (!this.progressBars[ds.name]){
                    return;
                }
                const n1 = ds.dataStore.size;
                const n2=  ds.dataStore.filterSize;
                this.progressBars[ds.name].setValue(n2);
                this.progressBars[ds.name].setText(n2);
                //this._callListeners(type,)
            }
        })
    }

    showAddViewDialog(){
        const controls =[
            {
                type:"checkbox",
                id:"clone-view",
                label:"Clone current view"
            },
            {
                type:"text",
                id:"name",
                label:"name"
            }

        ];
        if (this.dataSources.length>1){
            for (let ds of this.dataSources){
                controls.push({
                    type:"checkbox",
                    id:ds.name,
                    label:`Include ${ds.name}`


                })
                
            }
        }
        new CustomDialog({
            title:"Add New View",
            controls:controls,
            buttons:[{
                text:"Create New View",
                method:(vals)=>{
                    //create new view option
                    createEl("option",{text:vals["name"],value:vals["name"]},this.viewSelect);
                    this.viewSelect.value=vals["name"];
                    this.currentView=vals["name"];
                    if (!vals["clone-view"]){
                        //remove all charts and links
                        this.removeAllCharts();
                        this.viewData.links=[];
                        const state = this.getState();
                        const ic = state.view.initialCharts
                        for (let ds in this.dsIndex){
                            if (vals[ds]){
                                ic[ds]=[];
                            }
                        }
                        this._callListeners("state_saved",state);
                        this.contentDiv.innerHTML="";
                        this._init(state.view)
                    }
                    else{
                        const state = this.getState();
                        this._callListeners("state_saved",state);
                    }
                   
                    
                }
            }]
        })
    }

    showSaveViewDialog(action){
        new CustomDialog({
            title:"Save View",
            text:"Do you want to save the current view",
            buttons:[
            {
                text:"YES",
                method:()=>{
                    const state = this.getState();
                    this._callListeners("state_saved",state);
                    action();
                }
            },
            {
                text:"NO",
                method:()=>{
                    action();
                }
            }  
            ]
        })
    }


    changeView(view){
        this.removeAllCharts();
        this.contentDiv.innerHTML="";
        this.currentView=view;
        this.viewLoader(view).then(data=>{
            this._init(data);
        })
    }

    _columnRemoved(ds,col){      
        const ids_to_delete=[];
        for (let id in this.charts){
            const info = this.charts[id];
            if (info.dataSource===ds){
                const ch=info.chart;
                const div = ch.getDiv();
                const del = ch.onColumnRemoved(col);
                if (del){
                    div.remove(false);
                    ids_to_delete.push(id);
                    this._removeLinks(ch);
                    this._callListeners("chart_removed",ch);
                }
            }
        }
        //onColumnRemoved will remove the chart if it contains
        //data from the column, it will also remove the filter,
        //but not call any listeners
        if (ids_to_delete.length>0){
            ds.dataStore._callListeners("filtered"); 
        }
        for (let id of ids_to_delete){
            delete this.charts[id];
        }
    }

    _getColumnsRequiredForChart(config,set){
        const p = config.param;
      
        if (!p){
            return;
        }
        if (typeof p === "string"){
            set.add(p);
        }
        else{
            for (let i of p ){
                set.add(i);
            }
        }
        if (config.color_by){
            if (config.color_by.column){
                set.add(config.color_by.column.field);
            }
            else{
                set.add(config.color_by);
            }
            
        }
        if (config.tooltip){
            if (config.tooltip.column){
                set.add(config.tooltip.column);
            }
        }
        if (config.background_filter){
            set.add(config.background_filter.column);
        }
       
    }
   
    /**
    * Loads data from a remote file -the file must have headers (keys in the
    * case of json) which which match a columns field/id 
    * @param {object} info A config describing the file - 
    * @param {string} info.type - either csv,tsv ot json
    * @param {string} info.dataSource - the name of the datasource to load the data into
    * @param {string} info.url  - the url of the file
    * @param {function} [callback]  - a function to run once the data has loaded
    */
    loadFile(info,callback){
        const meths = {csv:csv,json:json,tsv:tsv}
        const iid  =  this.createInfoAlert("loading file",{spinner:true})
        meths[info.type](info.url).then(data=>{
            const cols={};
            const dataSource= info.dataSource;
            const ds =this.dsIndex[dataSource].dataStore;
            const all_cols =  ds.getAllColumns();
            //which columns are present in the datastore
            for (let c of data.columns){
                if (ds.columnIndex[c]){
                    cols[c]=[];
                }
            }
            for (let i=0;i<data.length;i++){
                const row =data[i];
                if (i+1%100===0){
                    this.updateInfoAlert(iid,`processed ${i}/${data.length} rows`);
                }
                for (let col in cols){
                    cols[col].push(row[col]);
                }
            }
            let proc=0;
            for  (let col in cols){
                ds.setColumnData(col,cols[col]);
                proc++;
                this.updateInfoAlert(iid,`processed ${proc}/${all_cols.length} columns`);       
            }
            this.updateInfoAlert(iid,"complete",{duration:2000})
            if (callback){
                callback();
            }
        })
    }


    _getUpdatedColumns(dataStore){
        const dc = dataStore.dirtyColumns;
        const rv = {
            columns:[],
            added:[],
            removed:[],
            colors_changed:[]
        }
        for (let c in dc.added){
            const td = getMd(c);
            rv.columns.push(td);
            rv.added.push(c)
        }
        for (let r in dc.removed){
            rv.removed.push(r)
        }

        for (let c in dc.data_changed){
            if (!rv.columns[c]){
                const td = getMd(c);
                rv.columns.push(td);
            }
        }

        for (let cc in dc.colors_changed){
            rv.colors_changed.push({
                column:cc,
                colors:dataStore.columnIndex[cc].colors
            })
        }

        return rv;
        
        function getMd(c){
            const cl = dataStore.columnIndex[c];
            const md={
                values:cl.values,
                datatype:cl.datatype,
                name:cl.name,
                editable:true,
                field:cl.field,

            }
           const arr = new Array(cl.data.length);
           for (let i=0;i<cl.data.length;i++){
               arr[i]= cl.data[i]
           }
           return {metadata:md,data:arr}

        }
    }


    getState(){
        const initialCharts={};
        const updatedColumns={};
        const metadata={};
        for (const ds of this.dataSources){
            if (ds.contentDiv){
                initialCharts[ds.name]=[];
            }
            updatedColumns[ds.name]=this._getUpdatedColumns(ds.dataStore); 
            const dstore= ds.dataStore;
            
            if (dstore.dirtyMetadata.size !==0){
                metadata[ds.name]={};
                for (let param of dstore.dirtyMetadata){
                    metadata[ds.name][param]=dstore[param];
                }
            }
        }
        for (let chid in this.charts){
            const chInfo = this.charts[chid];
           
            const chart = chInfo.chart;
            const config = chart.getConfig();
            const div =  chart.getDiv();
            config.position = [div.offsetLeft,div.offsetTop];
        
            initialCharts[chInfo.dataSource.name].push(config);
            
        }
       
        const view = JSON.parse(JSON.stringify(this.viewData))
        view.initialCharts= initialCharts;

        
        return{     
            view:view,
            currentView:this.currentView,
            updatedColumns:updatedColumns,
            metadata:metadata
        }
    }

    setAllColumnsClean(){
        for (let ds of this.dataSources){
            ds.dataStore.setAllColumnsClean();
        }
    }


    /** Displays a dialog
    * @param {Object} config extra settings
    * @param {Object[]} [config.initialCharts] A list of chart configs to
    * initially load
    * @param {function} [config.metaDataLoader] A method which returns a Promise
    * with a list of the column metadata, given the name of the datasource
    * and a list of column field/ids
    * @param {function} [config.onDataLoaded] A method which is called once all
    * the data is loaded (passing this as the only parameter)
    */

    showCustomDialog(config){
        new CustomDialog(config);
    }

     /**Adds a menu icon to either the main menubar or a datasource meubar
    * @param {string} datSource The name of dataa source or _main if addding
    * an icon to the main (top) toolbar
    * @param {string} icon The class name(s) of the icon
    * initially load
    * @param {string} text Text that will be diaplyed in a tooltip
    * @param {function} func The function that will be called when the icon is pressed
    */
    addMenuIcon(dataSource,icon,text,func){
        const pos = dataSource==="_main"?"bottom-right":"bottom";
        const el= dataSource==="_main"?this.leftMenuBar:this.dsIndex[dataSource].menuBar
        createMenuIcon(icon,{
            tooltip:{
                text:text,
                position:pos
            },
            func:func
        },el);
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


    /**
    * Loads data for specified columns into the appropriate dataStore
    * @param {string[]} columns An array of column fields/ids - if the datastore has no metadata
    * on a column, then it will call the metaDataLoader specified in the ChartManager's config
    * in order to load in the metadata, before loading in the actual data.
    * @param {string} dataSource The name of the dataSource
    * @param {function} callback A function which will be run once all the
    * columns are loaded
    * @param {integer} [split=10]  the number of columns to send with each request 
    * @param {integer} [threads=2]  the number of concurrent requests
    */
    loadColumnSet(columns,dataSource,callback,split=10,threads=2){

        //check to see if the datastore contains column information
        const ds = this.dsIndex[dataSource].dataStore;
        let noInfoCols = columns.filter(x=>!ds.columnIndex[x]);
        //load in the metadata
        if (noInfoCols.length>0){
            this.config.metaDataLoader(dataSource,noInfoCols).then((data)=>{
                for (let col of data){
                    ds.addColumn(col)
                }
                this._loadColumnSet(columns,dataSource,callback,split,threads);
            });
        }
        //metadata already present - load in the data 
        else{
            this._loadColumnSet(columns,dataSource,callback,split,threads);
        }
    }

    _loadColumnSet(columns,dataSource,callback,split=10,threads=2){

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
            this.columnsLoading[dataSource][col]=true;
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


    _loadColumnData(trans,dataSource){
        const dataStore=  this.dsIndex[dataSource].dataStore;
       
        const col_list = trans.columns[trans.nextColumn++];
        const columns=[];
        for (let col of col_list){
           columns.push(dataStore.getColumnInfo(col));
        }
        //float32 columns need to be at the beginning of the byte stream
        //as you can't create an array from  an arry buffer starting at
        //a byte position not divisible by 4 
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
            for (let col of col_list){
                delete this.columnsLoading[dataSource][col];
             }
            all_loaded = all_loaded>trans.totalColumns?trans.totalColumns:all_loaded;
            this.updateInfoAlert(trans.alertID,`Loading Columns:${all_loaded}/${trans.totalColumns}`);
            if (loaded>=total){
                this.updateInfoAlert(trans.alertID,`Loaded ${total} column${total===1?"":"s"}`,{duration:2000});
                trans.callback(trans.failedColumns);     
                delete this.transactions[trans.id];
            }
            if (trans.nextColumn<total){
                this._loadColumnData(trans,dataSource)
            }          
        })
    }


    _addLinkIcon(ds,ds_to,link){
        createMenuIcon("fas fa-plus-square",{
            tooltip:{
                text:`Add ${link.name}`,
                position:"bottom-right"
            },
            func:()=>{
                new AddColumnsFromRowsDialog(ds,ds_to,link,this);
            }
        },ds.menuBar);

    }


    _setUpMenu(ds){
        createMenuIcon("fas fa-chart-bar",{
            tooltip:{
                text:"Add Chart",
                position:"bottom-right"
            },
            func:()=>{
                new AddChartDialog(ds,config=>this.addChart(ds.name,config,true))
            }
        },ds.menuBar);

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
        if (ds.links){
            for (let ods in ds.links){
                const link= ds.links[ods];
                if (link.rows_as_columns){
                    this._addLinkIcon(ds,this.dsIndex[ods],link.rows_as_columns)
                }
            }

        }
        if (ds.link_to){
            createMenuIcon("fas fa-link",{
                tooltip:{
                    text:`Link (add) data to ${ds.link_to.dataSource} panel`,
                    position:"bottom-right"
                },
                func:()=>{
                   new LinkDataDialog(this,ds);
                }
    
            },ds.menuBar);
        }
        const idiv = createEl("div",{
            styles:{
                float:"right",
                lineHeight:"1.0"
            }
        },ds.menuBar);
        createEl("span",{
            text:ds.name,
            styles:{
                verticalAlign:"top",
                fontSize:"16px",
                marginRight:"4px"
            }
        },idiv);
        const size= ds.dataStore.size;
        ds.filterBar= createEl("progress",{
            value:size
        },ds.menBar)
        const pb = createEl("div",{
            styles:{
                width:"100px",
                display:"inline-block",
                marginTop:"2px"
            }
        },idiv);
        const pbConf={
            max:size,
            value:size,
            text:`${size}`
        }
        this.progressBars[ds.name]=new MDVProgress(pb,pbConf);
        
       
    }

    addListener(id,func){
        this.listeners[id]=func;
    }

    removeListener(id){
        delete this.listeners[id];
    }
    _callListeners(type,data){
        for (let id in this.listeners){
            this.listeners[id](type,this,data);
        }
    }

    /**
    * Adds a chart to the app
    * @param {string} dataSource The name of the chart's data source 
    * @param {string} config The chart's config
    * @param {boolean} [notify=false] If true any listeners will be informed that 
    * a chart has been loaded
    */
    addChart(dataSource,config,notify=false){
        //check if columns need loading
        const neededCols = new Set();
        this._getColumnsRequiredForChart(config,neededCols);
        //check which columns need loading
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
        const t = themes[this.theme];
        const div= createEl("div",{
            styles:{
                position:"absolute",
                width:width+"px",
                height:height+"px",
                left:left+"px",
                top:top+"px",
                background:t.main_panel_color,
                zIndex:2,
                display:"flex",
                alignItems:"center",
                justifyContent:"center"
            }
        },ds.contentDiv);
        createEl("i",{
            classes:["fas","fa-circle-notch","fa-spin"],
          
            styles:{
                fontSize:"30px",
                color:t.text_color
            }
        },div);
        createEl("div",{
            styles:{
                position:"absolute",
                overflow:"hide",
                textAlign:"center",
                top:"3px",
                color:t.text_color,
                textOverflow:"ellipsis",
                fontSize:"18px"

            },
            text:config.title
        },div)
        const func = ()=>{
            this._addChart(dataSource,config,div,notify);
        }
        this._getColumnsThen(dataSource,Array.from(neededCols),func);
    }

    
    _getColumnsFromOtherSource(dataSource,otherDataSource,columns,indexCol,func){
        this._getColumnsThen(otherDataSource,columns.concat(indexCol),()=>{
             const ds= this.dsIndex[dataSource].dataStore;
             const ods = this.dsIndex[otherDataSource].dataStore;
             const oindex = ods.getColumnIndex(indexCol);
             const ic = ds.columnIndex[indexCol]
             const index = ic.values.map(x=>oindex[x]);
             const colInfo= columns.map(x=>{
                const c1 = ds.columnIndex[x];
                const c2 = ods.columnIndex[x];
                if (c2.values){
                    c1.values=c2.values;
                }
                if (c2.minMax){
                    c1.minMax=c2.minMax;
                }
                if (c2.qauntiles){
                    c1.quantiles= c2.quantiles;
                }

                const buf = new  SharedArrayBuffer(ds.size * (c1.datatype==="text"?1:4));
                const arrType = c1.datatype==="text"?Uint8Array:Float32Array;
                return {
                    col:x,
                    data:buf,
                    arr:new arrType(buf),
                    odata:c2.data
                }
            });
            for (let n=0;n<ds.size;n++){
                const i = index[ic.data[n]];
                for (let c of colInfo){
                    c.arr[n]=c.odata[i]
                }
            }

            for (let c of colInfo){
                ds.setColumnData(c.col,c.data)
            }
            func();
        })

    }


    _getColumnsThen(dataSource,columns,func){
        const ds  =this.dsIndex[dataSource];
        const dStore = ds.dataStore;
        //check if need to load column data from linked data set
        if (ds.links){
            for (let ods in ds.links){
                const link = ds.links[ods];
                if (link.columns){
                    const otherCols = [];
                    const thisCols=[];
                    for (let c of columns){
                        if  (dStore.columnIndex[c].data){
                            continue;
                        }
                        if(link.columns.indexOf(c)===-1){
                            thisCols.push(c)
                        }
                        else{
                            otherCols.push(c)
                        }
                    }
                    //get the other datasource columns first
                    if (otherCols.length>0){
                        this._getColumnsFromOtherSource(dataSource,ods,
                            otherCols,link.index, ()=>{
                            this._getColumnsThen(dataSource,thisCols,func)
                        })
                        return;
                    }
                }
            }
        }
        const reqCols = columns.filter(x=>{
            //column already loading
            if (this.columnsLoading[dataSource][x]){
                return false;
            }
            const col = dStore.columnIndex[x];
            //no record of column- need to load it (plus metadata)
            if (!col){
                dStore.addColumnFromField(x);
                return true;
            }
            //only load if has no data
            return !col.data;
        });
      
        //No columns needed 
        //but columns requested by other actions may still be loading
        if (reqCols.length===0){
            this._haveColumnsLoaded(columns,dataSource,func);
        }
        //load requested columns then check all are loaded
        else{
            this.loadColumnSet(reqCols,dataSource,()=>{
                this._haveColumnsLoaded(columns,dataSource,func);
            })
        }
    }

    /*getIndexedData(dataSource,columns,indexColumn,callback,config={}){
        const col = this.dsIndex[dataSource].dataStore;
        this._getColumnsThen(dataSource,column,indexColumn],()=>{
            const index = ds.getColumnIndex(column);
            const cf = ds.getColorFunction(column,config);
            callback((val)=>{
                cf(index[val])
            })
        })

    }*/
    
    //need to ensure that column data is loaded before calling method
    _decorateColumnMethod(method,chart,dataSource){
        const newMethod = "_"+method;
        chart[newMethod]= chart[method];
        //if original method is called check whether column has data
        chart[method]=(column)=>{
            this._getColumnsThen(dataSource,[column],()=>chart[newMethod](column));
        }
    }

    //check all columns have loaded - if not recursive call after
    //time out, otherwise add the chart
    _haveColumnsLoaded(neededCols,dataSource,func){
        for (let col of neededCols){
            if (this.columnsLoading[dataSource][col]){
                setTimeout(()=>{
                    this._haveColumnsLoaded(neededCols,dataSource,func);
                },500);
                return;
            }
        }
        func();

    }

    _addChart(dataSource,config,div,notify=false){
        //**convert legacy data*********** 
        const ds= this.dsIndex[dataSource];
        div.innerHTML="";
        div.style.display="";
        div.style.alignItems="";
        div.style.justifyContent="";
        const chartType= BaseChart.types[config.type];
        const chart = new chartType.class(ds.dataStore,div,config);
        this._setChartTheme(chart);
        this.charts[chart.config.id]={
            chart:chart,
            dataSource:ds
        }
        this._makeChartRD(chart,ds);
        chart.popoutIcon = chart.addMenuIcon("fas fa-external-link-alt","popout",{
            func:()=>{
                this._popOutChart(chart,ds.contentDiv);
            }
        });     
        chart.addMenuIcon("fas fa-trash","remove chart")
            .addEventListener("click",()=>{   
                chart.remove();
                div.remove();
                delete this.charts[chart.config.id];
                this._removeLinks(chart);
                this._callListeners("chart_removed",chart);
            });
       
        //need to decorate any method that uses column data as data may
        //have to be loaded before method can execute
        if (chart.colorByColumn){
            this._decorateColumnMethod("colorByColumn",chart,dataSource);
        }
        if (chart.setToolTipColumn){
            this._decorateColumnMethod("setToolTipColumn",chart,dataSource);
        }
        if (chart.setBackgroundFilter){
            this._decorateColumnMethod("setBackgroundFilter",chart,dataSource);
        }
        if (chart.changeContourCategory){
            this._decorateColumnMethod("changeContourParameter",chart,dataSource);
        }

      

        const idl= ds.index_link_to;
        if (chart.setupLinks && idl){
            //ensures requested columns are loaded before other datasource loads them
            const func= (columns,callback)=>{
                //make sure index is loaded before use
                columns.push(idl.index);
                this._getColumnsThen(idl.dataSource,columns,callback)

            }    
            chart.setupLinks(this.dsIndex[idl.dataSource].dataStore, idl.index,func);
        }

        const cll= ds.column_link_to;
        if (cll && chart.createColumnLinks){
            const func= (columns,callback)=>{
                //make sure index is loaded before use
                this._getColumnsThen(cll.dataSource,columns,callback)
            }    
            chart.createColumnLinks(this.dsIndex[cll.dataSource].dataStore, cll.columns,func);
        }

        if (notify){
            this._callListeners("chart_added",chart);
        } 
        //check to see if all inital charted loaded , then can call any
        if (this._toLoadCharts){
            this._toLoadCharts.delete(config);
            if (this._toLoadCharts.size===0){
                delete this._toLoadCharts;
                if (this.viewData.links){
                    for (let l of this.viewData.links){
                        this._setUpLink(l);
                    }
                }
                this._callListeners("view_loaded",this.currentView)       
            }
        }
        return chart;
    }

    //sets up a link between charts
    _setUpLink(link){
        if (!link.id){
            link.id= getRandomString();
        }
        switch(link.type){
            case "color_by_column":
                const chart = this.charts[link.source_chart];
                if (!chart){
                    console.error(`broken link link:${link.id}`);
                }
                
              
                chart.chart.addListener(link.id,(type,data)=>{
                    if (type==="cell_clicked"){
                        for (let cid of link.target_charts){
                            this.getChart(cid).colorByColumn(data.row)
                        }
                    }
                })
                break;
        }
    }

    //if a chart has been removed, work out which links need to be removed
    _removeLinks(chart){
        const linksToRemove =[];
        const cid = chart.config.id;
        const links = this.viewData.links;
        if (!links){
            return;
        }
        for (let i=0;i<links.length;i++){
            const link= links[i];
            if (link.source_chart===cid){
                linksToRemove.push(i);
            }
            const index = link.target_charts.indexOf(cid)
            if (index!==-1){
                link.target_charts.splice(index,1);
                if (link.target_charts.length===0){
                    linksToRemove.push(i);
                }
            }
        }
        for (let i of linksToRemove){
            this.removeLink(i);
        }
    }

    removeLink(linkIndex){
        const link = this.viewData.links[linkIndex];
        switch(link.type){
            case "color_by_column":
                const chart =  this.charts[link.source_chart].chart;
                chart.removeListener(link.id)

        }
        this.viewData.links.splice(linkIndex,1)
    }

    

    removeAllCharts(dataSources){
        const allCharts=[];
        for (let cn in this.charts){
            const ch = this.charts[cn];
            if (dataSources && dataSources.indexOf(ch.dataSource.name) ===-1){
                continue;
            }
            allCharts.push([ch.chart,ch.window]);
        }
        for (let ci of allCharts){
            if (ci[1]){
                ci[1].close();

            }
            ci[0].remove()
            ci[0].div.remove()

        }
        this.charts={};
    }

    getAllFilters(dataSorce){
        const charts = this.getAllCharts(dataSorce);
        const fs= [];
        for (let c of charts){
            const filter = c.getFilter();
            if (filter){
                fs.push(filter)
            }
        }
        return fs;
    }

    

    getChart(id){
        const cinfo = this.charts[id];
        if (!cinfo){
            return null;
        }
        return cinfo.chart;
    }

    getAllCharts(dataSource){
        const charts= []
        for (let id in this.charts){
            const ch = this.charts[id];
            if (ch.dataSource.name ===dataSource){
                charts.push(ch.chart)
            }
        }
        return charts;
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
        },this.contentDiv)
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
                url:this.config.popouturl || "/"
        
            }
        );
        chInfo.window=win;
    }

    _sendAllChartsToBack(ds){
        for (let id in this.charts){  
            const c = this.charts[id]
            if (ds === c.dataSource ){
                c.chart.div.style.zIndex="";
            }
        }
    }

    _makeChartRD(chart,ds){
        const div = chart.getDiv();
        makeDraggable(div,{
            handle:".ciview-chart-title",
            contain:"topleft",
            ondragstart:(e)=>{
                this._sendAllChartsToBack(ds);
                div.style.zIndex=2;
            }
        });
        makeResizable(div,{
            onresizeend:(width,height)=>chart.setSize(width,height)
        })
    }
}



class SubGroupColumnChooser extends BaseDialog{
    constructor(dataStore,link,callback){
        const config={
            footer:true,
            width:250,
            maxHeight:500,
            title:"Select Columns",
            buttons:[{text:"OK",method:"getColumns"}]
        }
        super(config,{dataStore:dataStore,callback:callback,filter:filter});
        
    }
    init(content){
        createEl("input",{},this.div);
    }
}



/**
* Creates a dialog for the user to choose multiple columns
* @param {DataStore} dataStore - the dataStore the columns will be chosen from
* @param {function} callback - A function called when the user has selected the columns  
* The callback is provided with a list of chosen column fields(ids)
* @param {string} [filter=all] - The type of column the use can choose
*/

class ChooseColumnDialog extends BaseDialog{
    constructor(dataStore,callback,filter="all"){
        const config={
            footer:true,
            width:250,
            maxHeight:500,
            title:"Select Columns",
            buttons:[{text:"OK",method:"getColumns"}]
        }
        super(config,{dataStore:dataStore,callback:callback,filter:filter});
    }
    init(content){
        this.ds = content.dataStore;
        const gd= createEl("div",{styles:{padding:"8px"}});
        const rName = getRandomString();
        createEl("div",{text:"Groups"},this.dialog);
        
        const cgs = Object.keys(this.ds.columnGroups);
        cgs.unshift("All")
        for (let group of cgs){
            const d= createEl("span",{styles:{display:"inline-block",whiteSpace:"nowrap",marginRight:"5px"}},gd);
            createEl("span",{text:group},d)
            createEl("input",{
                type:"radio",
                value:group,
                name:rName
            },d)
            .addEventListener("click",e=>{
                this.checkAllInGroup(e.target.value);
            })
        }
        this.dialog.append(gd);
        createEl("div",{text:"Select Individual Columns"},this.dialog);
        const cd= createEl("div",{style:{padding:"8px"}});
        const cols = this.ds.getColumnList(content.filter);
        this.checks=[];
        this.callback=content.callback;
        for (let col of cols){
            const d= createEl("div",{
                styles:{//display:"inline-block",
                        whiteSpace:"nowrap",
                       // marginRight:"5px"
                    }
            },cd);
            //createEl("span",{text:col.name},d);
            const cb = createEl("input",{
                type:"checkbox"
            },d);
            this.checks.push([cb,col.field]);
            createEl("span",{text:col.name},d);
        }
        this.dialog.append(cd);
    }
    checkAllInGroup(group){
        if (group==="All"){
            for (let check of this.checks){ 
                check[0].checked=true;
            }
        }
        else{
            const cols = this.ds.columnGroups[group].columns;
            for (let check of this.checks){                  
                    check[0].checked=cols.indexOf(check[1])===-1?false:true
            }

        }
    
    }

    getColumns(){
        const cols=[];
        for (let check of this.checks){
            if (check[0].checked){
                cols.push(check[1]);
            }
        }
        this.callback(cols);
        this.close();
    }
}




/**
* Creates a dialog for the user to choose a chart and its associated parameters. When chosen the
* supplied callback will be invoked with the config of the chosen chart.
* @param {DataStore} dataStore - the dataStore the chart will be created from.
* @param {function} callback - A function called when the user has selected the chart and 
* its parameters. The callback is provided with the config of the chosen chart
*/
class AddChartDialog extends BaseDialog{
    constructor(dataSource,callback){
        const config={
            title:"Add Chart",
            columns:2,
            footer:true,
            width:380
        }
        super(config,{dataSource:dataSource,callback:callback});
        
    }
    init(content){
        this.extraControls={};
        const types=[];
        this.dataSource=content.dataSource;
        this.dataStore= content.dataSource.dataStore;
        for (let type in BaseChart.types){
            const t = BaseChart.types[type];
            //check to see if chart has any requirements
            let allow =true
            if (t.required){
                for (let r of t.required){
                    if (!this.dataSource[r]){
                        allow=false
                    }
                }
                if (!allow){
                    continue;
                }
            } 
            if (t.allow_user_add===false){
                continue;
            }
            types.push({
                name:t.name,
                type:type,
            });
        }
        
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
        const ed={};
        for (let name in this.extraControls){
            ed[name]= this.extraControls[name].value;
        }
        if (this.multiColumns){
            config.param =config.param.concat(this.multiColumns)
        }
        const t= BaseChart.types[this.chartType.value];

        if (t.init){
            t.init(config,this.dataSource,ed)
        }
        callback(config);
        this.chartName.value="";
        this.chartDescription.value="";
        this.chartType.value= this.defaultType;
        this.setParamDiv(this.defaultType)

    }

    _addMultiColumnSelect(holder,filter){
        //get default values
        const ps= this.dataStore.getColumnList(filter);
        let text ="";
        if (ps.length>1){
            text= `${ps[0].name},... (1)`
            this.multiColumns=[ps[0].field];
        }
        const dd = createEl("span",{text:text},holder);
        createEl("i",{classes:["fas","fa-plus"]},holder)
        holder.style.cursor="pointer";
        holder.addEventListener("click",()=>{
            new ChooseColumnDialog(this.dataStore,cols=>{
                this.multiColumns=cols;
                let text="";
                if (cols.length>0){
                    const max= cols.length<3?cols.length:3;
                    const arr= []
                    for (let n=0;n<max;n++){
                        arr.push(this.dataStore.getColumnName(cols[n]));
                    }
                    text = arr.join(",");
                    if (cols.length>3){
                        text+=",...."
                    }
                    text+=`(${cols.length})`;
                    dd.textContent=text;
                }
               
            },filter);
        });
    }

    setParamDiv(type){
        this.paramDiv.innerHTML="";
        const params = BaseChart.types[type].params;
        this.paramSelects=[];
        if (params){
            for (let p of params){
                const d = createEl("div",{styles:{padding:"4px"}},this.paramDiv)
                const sp =createEl("div",{text:p.name+":"},d);
                const holder =createEl("div",{},this.paramDiv);
                if (p.type.startsWith("_multi")){
                this._addMultiColumnSelect(holder,p.type.split(":")[1])
                }
                else{
                    this.multiColumns=null;
                    const dd = createEl("select",{
                        styles:{
                            maxWidth:"200px"
                        }
                    },holder);
                    const ps= this.dataStore.getColumnList(p.type);
                    const sgs = {}
                    for (let ds of this.dataStore.subgroupDataSources){
                        sgs[ds]=createEl("optgroup",{label:ds});
                    }
                    for (let item of ps){
                        let ele = dd;
                        if (item.subgroup){
                            ele=sgs[item.subgroup.dataSource];
                        }
                        createEl("option",{text:item.name,value:item.field},ele);

                    }
                    for (let ds of this.dataStore.subgroupDataSources){
                        dd.append(sgs[ds]);
                    }
                    this.paramSelects.push(dd);
                }           
            }
        }
        const t= BaseChart.types[this.chartType.value];
        this.extraControls={};
        if (t.extra_controls){
            const controls = t.extra_controls(this.dataSource);
            for (let c of controls){
                createEl("div",{
                    text:c.label,
                    classes:["ciview-title-div"]
                },this.columns[1]);
                if (c.type==="dropdown"){
                    const sel = createEl("select",{
                        styles:{
                            maxWidth:"200px"
                        }
                    },this.columns[1]);
                    
                    for (let item of c.values){
                        createEl("option",{text:item.name,value:item.value},sel)
                    }
                    this.extraControls[c.name]=sel;
                }
            }

        }
    }
}

export default ChartManager;
export {AddChartDialog,ChooseColumnDialog};