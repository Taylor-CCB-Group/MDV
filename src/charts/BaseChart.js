import {getRandomString} from "../utilities/Utilities.js";
import {ContextMenu} from "../utilities/ContextMenu.js";
import {createEl} from "../utilities/Elements.js";
import {ChartSettingsDialog} from "./ChartSettingsDialog.js";


class BaseChart{
    /**
     * The base constructor
     * @param {DataStore} dataStore - The datastore object thta contains the data for this chart
     * @param {string} div - The id of the div element or the element itself to house the graph
     * @param {Object} config - The config describing the graph
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
        });
   
        //set up context menu and icon which opens it
        this.contextMenu = new ContextMenu((data)=>{     
            return this.getContextMenu(data);
        })

        this.addMenuIcon("fas fa-bars","more")
        .addEventListener("click",e=>this.contextMenu.show(e));
        
        this.addMenuIcon("fas fa-cog","settings")
        .addEventListener("click",(e)=>this._openSettingsDialog(e));

      
           
        //info icon
        this.addMenuIcon("fas fa-info",config.legend || "No description",{size:"medium"});
     
        this.div.append(this.titleBar);
        this.div.append(this.contentDiv); 
        this.listeners={};

        //work out width and height based on container
        this._setDimensions();
    }


    /**
    * Adds a listener to the datastore that will be called when an event occurs,
    * passing the event type and any data. There are the following different types
    * of event:-
    * <ul>
    * <li> removed - called when thr chart has been removed </li>
    * </ul>
    * @param {string} id- a unique id indetifying the listener
    * @param {function} listener - a function that accepts two paramaters, the type
    * of event and the dat associated with it
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
    * @param {string} icon- the css classs of the icon (space delimited).
    * @param {string} tootltip - the tooltip text
    * @param {object} config - extra inormation about the icon/tooltip
    * @param {string} [config.size=small] - the size of the tooltip
    * @param {string} [config.position=bottom] - the position of the tooltip.
    * @param {function} [config.func=] - a function that is called when the icon is clicked
    * @returns {DOMElement} - the icon
    */
    addMenuIcon(icon,tooltip,config={}){
        const sp= createEl("span",{
            "aria-label":tooltip,
            role:"tooltip",
            "data-microtip-size":config.size || "small",
            "data-microtip-position":config.position || "bottom-left",
            styles:{
                margin:"0px 1px"
            }
            },this.menuSpace);
    
        createEl("i",{  
            classes:["ciview-chart-icon"].concat(icon.split(" "))
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
    * be implemented on any subclasses
    */
    onDataFiltered(){}

 

    /**
    * Just removes the DOM elements, subclasses should do their own cleanup
    */
    remove(){
        this.titleBar.remove()
        this.contentDiv.remove();
        this.dataStore.removeListener(this.config.id);
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
            const cols = this.dataStore.getColumnList(colorOptions.colorby);
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
            settings.push({
                label:"Show Color Legend",
                type:"check",
                values:[cols,"name","field"],
                current_value:c.color_legend?c.color_legend.display:true, 
                func:(x)=>{
                    if (!c.color_by){
                        return
                    }
                    c.color_legend.display=x;
                    this.addColorLegend();

                }

            })
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
        if (!this.settingsDialog){
          this.settingsDialog=  new ChartSettingsDialog({
                maxHeight:400,
                doc:this.__doc__ || document,
                width:300,
                title:"Settings",
                position:[e.pageX,e.pageY],
                onclose:()=>this.settingsDialog=null

            },this.getSettings());

        }
      
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
        if (this.settingsDialog){
            this.settingsDialog.close();
        }
        if (this.legend){
            this.legend.__doc__=doc;
        }
    }

    _setConfigValue(conf,value,def){
        if (conf[value]===undefined){
            conf[value]=def;
        }
        return conf[value]
    }



    
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
        this._setDimensions()
    }

    _setDimensions(){
        const rect = this.div.getBoundingClientRect();
        const y=rect.height;
        const x=rect.width;
        this.config.size=[x,y];
        this.height=y-20;
        this.width=x;
    }


    /**
     * Downloads an image of the chart
     * @param {string} im_type - either svg or png 
     */
    downloadImage(im_type){    
        this.getImage(resp=>{
            let link =document.createElement("a");
            link.download=this.config.title+"."+im_type;
            if (im_type==="svg"){
                link.href=resp;
            }
            else{
                let url =resp.toDataURL('image/png');
                url = url.replace(/^data:image\/png/,'data:application/octet-stream');
                link.href=url
            }        
            link.click();
            link.remove()
        },im_type);       
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

BaseChart.types={};

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
 

  

export {BaseChart}
