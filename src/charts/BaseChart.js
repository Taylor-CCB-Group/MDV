import {getRandomString} from "../utilities/Utilities.js";
import {ContextMenu} from "../utilities/ContextMenu.js";
import {createEl} from "../utilities/Elements.js";
import SettingsDialog from "../utilities/SettingsDialog";
import { chartTypes } from "./ChartTypes";
import DebugJsonDialogReactWrapper from "../react/components/DebugJsonDialogReactWrapper";


class BaseChart{
    /**
     * The base constructor
     * @param {DataStore} dataStore - The datastore object that contains the data for this chart
     * @param {string | HTMLDivElement} div - The id of the div element or the element itself to house the chart
     * @param {Object} config - The config describing the chart
     */
    constructor(dataStore,div,config){
        //******adapt legacy configs
        if (config.color_by){
            if (config.color_by.column){
                config.color_by = config.color_by.column.field;
            }
        }
        //**********


        //copy the config
        this.config=JSON.parse(JSON.stringify(config))
        //give it a random id if one isn't supplied
        if (!this.config.id){
            this.config.id=getRandomString();
        }
        //required in case added to separate browser window
        this.__doc__ = document;

        this.dataStore=dataStore;
        this.listeners={};
        
        //create the DOM elements
        this.div = typeof div === "string"?document.getElementById(div):div;
        this.div.classList.add("ciview-chart-panel");
        this.titleBar=createEl("div",{
            classes:["ciview-chart-titlebar"]
        })
        
        this.title= createEl("div",{
            classes:["ciview-chart-title"],
            text:config.title,
        },this.titleBar)

        this.menuSpace = createEl("div",{
            classes:["ciview-chart-menuspace"]
        },this.titleBar);
        
        this.contentDiv = createEl("div",{
            classes:["ciview-chart-content"]
        });

        if (config.title_color){
            this.titleBar.style.backgroundColor=config.title_color;
        }

        //reset button
        this.resetButton = this.addMenuIcon("fas fa-sync","remove filter");
        this.resetButton.style.display="none";
        this.resetButton.addEventListener("click",()=>{
            this.removeFilter();
            this.resetButton.style.display="none";
        });



      
        //register with datastore to listen to filter events
        this.dataStore.addListener(this.config.id,(type,data)=>{
            if (type==="filtered"){
                this.onDataFiltered(data)
            }
            else if (type==="data_changed"){
                this.onDataChanged(data);
            }
            else if (type==="data_added"){
                this.onDataAdded(data);
            }
            else if (type==="data_highlighted"){
                if (data.source === this) {
                    this._callListeners("data_highlighted", data);
                }
                if (this.onDataHighlighted){
                    this.onDataHighlighted(data);
                }
            }

        });
   
        //set up context menu and icon which opens it
        this.contextMenu = new ContextMenu((data)=>{     
            const menu = this.getContextMenu(data);
            if (import.meta.env.MODE !== "production") {
                menu.push({
                    text:"debug chart",
                    icon:"fas fa-bug",
                    func:()=>{
                        window.mdv.debugChart = this;
                        this.dialogs.push(new DebugJsonDialogReactWrapper(this.config, this));
                    }
                });
                menu.push({
                    text: "copy config JSON to clipboard",
                    icon: "fas fa-copy",
                    func: () => navigator.clipboard.writeText(JSON.stringify(this.config, null, 2))
                });
            }
            return menu;
        })

        this.addMenuIcon("fas fa-bars","more")
        .addEventListener("click",e=>this.contextMenu.show(e));
        
        this.addMenuIcon("fas fa-cog","settings")
        .addEventListener("click",(e)=>this._openSettingsDialog(e));

      
           
        //info icon
        this.legendIcon = this.addMenuIcon("fas fa-info",config.legend || "No description",{size:"medium"});

        let oldSize = config.size;
        this.contentDiv.addEventListener("fullscreenchange", ()=>{
            //nb, debounced version of setSize also being called by gridstack - doesn't seem to cause any problems
            if (document.fullscreenElement) {
                if (this.contentDiv !== document.fullscreenElement) console.error('unexpected fullscreen element');
                const rect = this.contentDiv.getBoundingClientRect();
                this.setSize(rect.width, rect.height);
                for (const d of this.dialogs) {
                    d.setParent(this.contentDiv);
                }
            } else {
                this.setSize(...oldSize);
                for (const d of this.dialogs) {
                    d.setParent(null);
                }
            }
        });
        this.addMenuIcon("fas fa-expand","fullscreen", {
            func: async ()=>{
                oldSize = this.config.size;
                await this.contentDiv.requestFullscreen();
            }
        });
     
        this.div.append(this.titleBar);
        this.div.append(this.contentDiv); 
        this.listeners={};

        //work out width and height based on container
        this._setDimensions();
    }
    dialogs = [];
    _getContentDimensions(){
        return{ //PJT to review re. gridstack.
            top:5,
            left:5,
            height:this.height,
            width:this.width-5
        }
    }

    getFilter(){

    }

    setFilter(){
        
    }


    /**
    * Adds a listener to the datastore that will be called when an event occurs,
    * passing the event type and any data. There are the following different types
    * of event:-
    * <ul>
    * <li> removed - called when the chart has been removed </li>
    * </ul>
    * @param {string} id - a unique id indetifying the listener
    * @param {function} listener - a function that accepts two paramaters, the type
    * of event and the data associated with it
    */
    addListener(id,listener){
        this.listeners[id]=listener;
    }

    /**
    * Removes the specified listener from the chart
    * @param {string} id The id of the listener to remove 
    */
    removeListener(id){
        delete this.listeners[id];
    }

    _callListeners(type,data){
        for (let id in this.listeners){
            this.listeners[id](type,data);
        }
    }

    /**
    * Adds a menu icon with tooltip to the title bar 
    * @param {string} icon - the css classs of the icon (space delimited).
    * @param {string} tooltip - the tooltip text
    * @param {object} config - extra inormation about the icon/tooltip
    * @param {string} [config.size=small] - the size of the tooltip
    * @param {string} [config.position=bottom] - the position of the tooltip.
    * @param {function} [config.func] - a function that is called when the icon is clicked
    * @returns {DOMElement} - the icon
    */
    addMenuIcon(icon,tooltip,config={}){
        const sp = createEl("span", {
            "aria-label": tooltip,
            "data-microtip-color": "red",
            role: "tooltip", 
            "data-microtip-size": config.size || "small",
            "data-microtip-position": config.position || "bottom-left",
            styles: {
                margin:"0px 1px"
            }
        }, this.menuSpace);
    
        createEl("i",{  
            classes:["ciview-chart-icon"].concat(icon.split(" ")) //a11y - we could use an actual button
        },sp);
        if (config.func){
            sp.addEventListener("click",(e)=>config.func(e));
        }
        return sp;

    }

    /**
    * needs to be implemented by subclasses
    */
    removeFilter(){}

    /**
    * Called by the datastore when the data is filtered. Needs to
    * be implemented on any subclasses.
    * @param {Dimension} dim - the dimension that has been filtered
    */
    onDataFiltered(dim){}

    /**Check if chart is composed of any columns whose data has
     * changed. if so re-calculate and re-draw chart (call onDataFiltered)
     * @param {object} data - a list of column/fields whose data has been modified
     * @param {string[]} data.columns a list of column ids whose data has changed
     * @param {boolean} data.hasFiltered Whether a 'filtered' callback has already been
     * issued 
     */
    onDataChanged(data){
        const columns = data.columns;
        //update any charts which use data from the columns
        //(if they haven't already been updated by the filter changing)
        if (!data.hasFiltered){
            let cols= this.config.param;
            let isDirty=false
            if (typeof this.config.param === "string" ){
                cols=[this.config.param];
            }
            for (let p of cols){
                if (columns.indexOf(p)!==-1){
                    isDirty=true;
                    break;
                }
            }
            if (isDirty){
                this.onDataFiltered();
            }

        }
        //recolor any charts coloured by the column
        if (columns.indexOf(this.config.color_by)!==-1){
            this.colorByColumn(this.config.color_by);
        }
    }





    getColorLegend(){
        const conf = {overideValues:{
            colorLogScale:this.config.log_color_scale
        }};
        this._addTrimmedColor(this.config.color_by,conf);

        return this.dataStore.getColorLegend(this.config.color_by,conf);
    }


    getQunatile

    _addTrimmedColor(column,conf){
        const tr = this.config.trim_color_scale;
        const col= this.dataStore.columnIndex[column];
        if (tr && tr !=="none"){
            if (col.quantiles && col.quantiles !=="NA"){
                conf.overideValues.min=col.quantiles[tr][0];
                conf.overideValues.max=col.quantiles[tr][1];
            }
        }

    }

    /**
     * adds (or removes) the color legend depending on the chart's
     * config color_legend.display value - assumes chart has a 
     * get colorLegend method 
     */
     setColorLegend(){
        if (!this.config.color_legend.display){
            if (this.legend){
                this.config.color_legend.pos=[this.legend.offsetLeft,this.legend.offsetTop];
                this.legend.remove();
                delete this.legend;
            }
            return;
        }
        const box = this._getContentDimensions();
        let lt= 0;
        let ll = 0;
        if (this.legend){
            ll = this.legend.style.left;
            lt = this.legend.style.top;
            this.legend.remove();
        }
        else {
            const cl= this.config.color_legend;
            if (!cl.pos){
                cl.pos=[box.left,box.top]
            }
            ll= cl.pos[0]+"px";
            lt= cl.pos[1]+"px";
        }
        this.legend = this.getColorLegend();
        if (!this.legend) {
            console.warn('no color legend');
            return;
        }
        this.contentDiv.append(this.legend);
       
        this.legend.style.left= ll;
        this.legend.style.top= lt;
        this.legend.__doc__=this.__doc__;

    }

    getColorFunction(column,asArray){
        this.config.color_by=column;
        const conf ={
            asArray:asArray,
            overideValues:{
                colorLogScale:this.config.log_color_scale,
                fallbackOnZero: this.config.fallbackOnZero
            }
        }
        this._addTrimmedColor(column,conf);
       
      
		const colorFunc = this.dataStore.getColorFunction(column,conf);
       
        if (!this.config.color_legend){
            this.config.color_legend={
                display:true
            }
        }
        
        this.setColorLegend();
        return colorFunc;
        
    }


    /**Checks to see if the column is used in the chart
    * If so, the chart will be removed but no callbacks will be involved
    * @param {object} data - a list of column/fields whose data has been modified
    * @param {string[]} data.columns a list of column ids whose data has changed
    * @param {boolean} hasFiltered Whether a 'filtered' callback has already been
    * issued 
    * @returns {boolean} true if the chart has been removed
    */
    onColumnRemoved(column){
        let cols= this.config.param;
        let isDirty=false
        if (typeof this.config.param === "string" ){
            cols=[this.config.param];
        }
        for (let p of cols){
            if (column===p){
                isDirty=true;
                break;
            }
        }
        if (isDirty){
            this.remove(false);
            return true;
        }
        if (this.colorByColumn){
            if (this.config.color_by===column){
                delete this.config.color_by;
                this.colorByDefault();
            }
        }
        return false;

    }

    onDataAdded(newSize){
        this.onDataFiltered();
    }

    onDataHighlighted(data){}

    addToolTip(){
        this._tooltip = createEl("div",{
            classes:["ciview-tooltip"],
            stlyles:{
                display:"none",
                position:"absolute"
            }
        },this.__doc__.body)
    }

    showToolTip(e,msg){
        if (!this._tooltip) this.addToolTip();
        this._tooltip.innerHTML=msg;
        this._tooltip.style.display= "inline-block";
        this._tooltip.style.left= (3+e.clientX)+"px";
        this._tooltip.style.top=(3+e.clientY)+"px"
    }

    hideToolTip(){
        if (!this._tooltip) return;
        this._tooltip.style.display="none";
    }

    /**
    * Just removes the DOM elements, subclasses should do their own cleanup
    */
    remove(){
        this.titleBar.remove()
        this.contentDiv.remove();
        this.dataStore.removeListener(this.config.id);
        if (this._tooltip){
            this._tooltip.remove();
        }
        for (const d of this.dialogs){
            d.close();
        }
        // dynamic props?
    }


     /**
    * Returns information about which controls to add to the settings dialog.
    * Subclasses should call this method and then add their own controls e.g.
    * <pre style='background:lightgray;padding:10px'>
    * getSettings(){
    *     let settings = super.getSettings();
    *     return settings.concat([{
    *       label:"A value"
    *       type:"slider",
    *       default_value:this.config.value,
    *       max:10,
    *       min:10,
    *       func:x=>{
    *           this.config.value=x;
    *           //update chart
    *       }      
    *     }]);
    * }
    * </pre>
    */
    getSettings(){
        const c= this.config;
        let settings = [
            {
                type:"text",
                label:"Chart Name",
                current_value:c.title,
                func:(v)=>{
                    this.setTitle(v);
                }
            }
        ];
        const colorOptions = this.getColorOptions();
        
        if (colorOptions.colorby){
            //cannot color by unique (at the moment)
            const filter = colorOptions.colorby==="all"?["int32","text","integer","double","text16","multitext"]:colorOptions.colorby;
            const cols = this.dataStore.getColumnList(filter);
            cols.push({name:"None",field:"_none"})
            settings.push({
                label:"Color By",
                type:"dropdown",
                values:[cols,"name","field"],
                current_value:c.color_by || "_none",
                func:(x)=>{
                    if (x==="_none"){
                        delete c.color_by
                        this.colorByDefault();
                    }
                    else{
                        c.color_by=x;
                        this.colorByColumn(x);
                    }

                }
            });
            if (colorOptions.color_overlay !== undefined) {
                settings.push({
                    label:"Color Overlay",
                    type:"slider",
                    current_value:c.color_overlay,
                    func:(x)=>{
                        c.color_overlay=x;
                        this.colorByColumn(c.color_by);
                    }
                });
            }
            settings.push({
                label:"Show Color Legend",
                type:"check",
            
                current_value:c.color_legend?c.color_legend.display:true, 
                func:(x)=>{
                    if (!c.color_by){
                        return
                    }
                    c.color_legend.display=x;
                    this.setColorLegend();

                }

            });
            settings.push({
                label:"SymLog Color Scale",
                type:"check",
                
                current_value:c.log_color_scale, 
                func:(x)=>{      
                    c.log_color_scale=x;
                    if (c.color_by){
                        this.colorByColumn(c.color_by);
                    }
                }
            });
            settings.push({
                label:"Treat zero as missing",
                type:"check",
                
                current_value:c.fallbackOnZero,
                func:(x)=>{      
                    c.fallbackOnZero = x;
                    if (c.color_by){
                        this.colorByColumn(c.color_by);
                    }
                }
            });
            settings.push({
                type:"radiobuttons",
                label:"Trim Color Scale to Percentile",
                current_value:c.trim_color_scale || "none",
                choices:[["No Trim","none"],["0.001","0.001"],["0.01","0.01"],["0.05","0.05"]],             
                func:(v)=>{
                    c.trim_color_scale=v;
                    if (c.color_by){
                        this.colorByColumn(c.color_by);
                    }
                   
                }
            });

        }

        return settings;
             
    }



    _addUnpinIcon(){
        this.unpinIcon = this.addMenuIcon("fas fa-thumbtack","unpin chart")
        this.unpinIcon.addEventListener("click",(e)=>{
            this.unpinIcon.remove();
            this.unpinChart();

        })
    }

   
    _openSettingsDialog(e){
        if (!this.settingsDialog) {
            if (import.meta.env.DEV) {
                // this.settingsDialog = new SettingsDialogReactWrapper(this);
                // return;
            };
            this.settingsDialog = new SettingsDialog({
                maxHeight:400,
                doc:this.__doc__ || document,
                width:300,
                title:"Settings",
                position:[e.pageX,e.pageY],
                useMobx: this.useMobx,
                onclose:()=>this.settingsDialog=null
            },this.getSettings());
            this.dialogs.push(this.settingsDialog);
        }
        //experimenting with making the parent always be contentDiv (not only in fullscreen mode)
        //doesn't work ATM because of stacking context - may want to review that more generally 
        //(can be annoying with tooltips...)
        //this.settingsDialog.setParent(this.contentDiv);
    }

  
    getContextMenu(data){
        let menu= [];
        if (this.getImage){
            menu=menu.concat([ 
            {
               text:"create svg image",
               icon:"far fa-image",
               func:()=>this.downloadImage("svg")
            },
            {
                text:"create png image",
                icon:"far fa-image",
                func:()=>this.downloadImage("png")
            }
            ]);
        }          
        if (this.pinChart && !(this.isPinned)){
            menu.push({
                text:"pin chart",
                icon:"fas fa-thumbtack",
                func:()=>{
                    this.pinChart();
                    this._addUnpinIcon();
                }
            })
        }

        if (this.getChartData){
            menu.push({
                text:"download data",
                icon:"fas fa-download",
                func:()=>this.downloadData()
            })
        }
        if (this.addToContextMenu){
            menu=menu.concat(this.addToContextMenu());
        }
        menu.push({
            text: "experimental settings dialog",
            icon: "fas fa-cog",
            func: async () => {
                const m = await import("../react/components/SettingsDialogReactWrapper");
                const SettingsDialogReactWrapper = m.default;
                this.dialogs.push(new SettingsDialogReactWrapper(this));
            }
        });
        return menu;
    }


    /**
    * Returns the dom element that the chart is attached to
    */
    getDiv(){
        return this.div;
    }

    /**
    * Instructs the chart to use a different document. This is only required if you are 
    * going to add the chart to a different browser window
    * @param {document} doc - the document that the chart will use 
    */
    changeBaseDocument(doc){
        this.contextMenu.__doc__=doc;
        this.__doc__=doc;
        for (const d of this.dialogs){
            d.close();
        }
        if (this.legend){
            this.legend.__doc__=doc;
        }
        if (this._tooltip){
            this._tooltip.remove();
            this.addToolTip();
        }
        if (this.extra_legends){
            for (let l of this.extra_legends){
                if (this[l]){
                    this[l].__doc__=doc;
                }
            }
        }
    }

    _setConfigValue(conf,value,def){
        if (conf[value]===undefined){
            conf[value]=def;
        }
        return conf[value]
    }

    drawChart(){}



    
    /**
    * Returns a copy of the chart's config
    */
    getConfig(){
        if (this.legend){
            this.config.color_legend={
                display:true,
                pos:[this.legend.offsetLeft,this.legend.offsetTop]
            }
        }
        return JSON.parse(JSON.stringify(this.config))
    }

    getColorOptions(){
        return {};
    }


    /**
     * Sets the size of the graph. If no parameters are supplied
     * then the graph will be resized based on it container. 
     * Subclasses should overide this, but call the super method
     * first which will calculate width and height of the content div
     * @param {integer} x - The new width 
     * @param {integer} y The new height;
     */
    setSize(x,y){
        //if supplied change the div dimensions
        if (x){
            this.div.style.height=y+"px";
            this.div.style.width=x+"px";
        }
        //calculate width and height based on outer div
        this._setDimensions();
    }

    _setDimensions(){
        const rect = this.div.getBoundingClientRect();
        const y=Math.round(rect.height);
        const x=Math.round(rect.width);
        this.config.size=[x,y];
        const rect2  = this.contentDiv.getBoundingClientRect();
        this.height=Math.round(rect2.height);
        this.width=x;
    }


    /**
     * Downloads an image of the chart
     * @param {string} im_type - either svg or png 
     */
    downloadImage(im_type){ 
        const originalColor =this.contentDiv.style.color;
        this.contentDiv.style.color = "black";
        this.getImage(resp=>{
            let link =document.createElement("a");
            let name = this.config.title || "image"
            link.download=name+"."+im_type;
            if (im_type==="svg"){
                link.href="data:image/svg+xml," + encodeURIComponent(resp);
            }
            else{
                let url =resp.toDataURL('image/png');
                url = url.replace(/^data:image\/png/,'data:application/octet-stream');
                link.href=url
            }        
            link.click();
            link.remove();
            this.contentDiv.style.color = originalColor;
        },im_type);       
    }

    downloadData(){
        const blob = this.getChartData();
        const save = createEl("a",{
            download:this.config.title,
            target:"_blank",
            href:window.URL.createObjectURL(blob)

        },this.__doc__.body)                
        save.click();
        save.remove();
    }

  
  setTitle(title){
       this.title.textContent=title;
       this.config.title=title;

  }

  getImageFromSVG(svg,callback) {
    var copy = svg.cloneNode(true);
    copyStylesInline(copy, svg);
    var canvas = document.createElement("canvas");
    //var bbox = svg.getBBox();
    copy.style.top = "0px";
    canvas.width = svg.width.baseVal.value
    canvas.height =svg.height.baseVal.value
    var ctx = canvas.getContext("2d");
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    var data = (new XMLSerializer()).serializeToString(copy);
    var DOMURL = window.URL || window.webkitURL || window;
    var img = new Image();
    var svgBlob = new Blob([data], {type: "image/svg+xml;charset=utf-8"});
    var url = DOMURL.createObjectURL(svgBlob);
    img.src = url;
          img.onload = function () {
          ctx.drawImage(img, 0, 0);
          callback(canvas,ctx)
    }
   
  }


  


}

/**
 * A dictionary of all the chart types
 */
BaseChart.types = chartTypes;

function copyStylesInline(destinationNode, sourceNode) {
    var containerElements = ["svg","g"];
    for (var cd = 0; cd < destinationNode.childNodes.length; cd++) {
        var child = destinationNode.childNodes[cd];
        if (containerElements.indexOf(child.tagName) != -1) {
             copyStylesInline(child, sourceNode.childNodes[cd]);
             continue;
        }
        var style = sourceNode.childNodes[cd].currentStyle || window.getComputedStyle(sourceNode.childNodes[cd]);
        if (style == "undefined" || style == null) continue;
        for (var st = 0; st < style.length; st++){
             child.style.setProperty(style[st], style.getPropertyValue(style[st]));
        }
    }
}
 

  

export default BaseChart
