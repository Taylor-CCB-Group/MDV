import BaseChart from "./BaseChart.js";
import {createEl} from "../utilities/Elements.js";
import CustomDialog from "./dialogs/CustomDialog.js";




class GenomeBrowser extends BaseChart{
    constructor(dataStore,div,config){
		super(dataStore,div,config);
        const c = this.config;
        c.view_margins= c.view_margins || {type:"percentage",value:20};
        this.contentDiv.style.width="100%";
        this.contentDiv.style.height="100%";
        c.type="genome_browser";
        const add_ruler= c.tracks.find(x=>x["format"]==="ruler")?false:true;

        import ('../browser/panel.js?v0.1').then(({default:MLVPanel})=>{
            this.browser = new MLVPanel(this.contentDiv,{
                allow_user_interactions:true,
                show_scale:true,
                ruler_track:add_ruler
            },c.tracks
            );
    
            this.bamscatrack=this.browser.getTrack("_atac_bam_track");
            this.baseTrack= this.browser.getTrack("_base_track");
            this.baseTrack.dataStore = this.dataStore;
       
            this.addMenuIcon("fa fa-eye-slash",config.legend || "Shoow/Hide Tracks",
                {
                    func:e=>this.showHideDialog()
    
                });
            this.locText=createEl("input",{styles:{fontSize:"11px"}},this.menuSpace);
            this.locText.addEventListener("keypress",e=>{
                if (e.keyCode==13){
                    const loc = this._calculatePosition(this.locText.value);
                    this.browser.update(loc.chr,loc.start,loc.end);
                }
            })
            this.browser.addListener(this.config.id,(type,data)=>this.onBrowserAction(type,data));
            if (c.feature_label){
                this.setLabelFunction(c.feature_label);
            }
            if (c.color_by){
                 this.colorByColumn(c.color_by,false)
            }
            if (c.color_wig_tracks_by){
                const col = c.color_wig_tracks_by;
                const vc= this.dataStore.getValueToColor(col);
                //color any wig tracks
                for (let v in vc){
                    const tr =this.browser.tracks[`${col}|${v}`];
                    if (tr){
                        tr.config.color= vc[v];
                    }
                }
            }    
            if (!this.bamscatrack){
                const g= this.config.genome_location;
                if (g){
                    
                    this.browser.update(g.chr,Math.round(g.start),Math.round(g.end),true);
                }
                else{
                    this.onDataHighlighted({indexes:[0]});
                }
            }
        });
	}



    onBrowserAction(type,data){
        switch(type){
            case "featureclick":
                if (data.track.config.track_id==="_base_track"){
                    const fIndex = data.feature.data[0];
                    this.dataStore.dataHighlighted([parseInt(fIndex)],this);
                }
                if (data.track.config.track_id==="_bam"){
                    console.log(data.feature);                
                }
                break;
            case "range_selected":
                if (this.bamscatrack){
                    const ids =this.bamscatrack.getIdsInRange(data);
                    this.cellDim.filter("filterOnIndex",[],ids);
                    this.browser.setHighlightedRegion(data,"_filter","blue");
                    this.resetButton.style.display="inline";
                    this.browser.update();         
                }
                break;
            case "view_changed":
                this.locText.value=`${data.chr}:${data.start}-${data.end}`;

        }
    }

    //called when chart is created passing the linked datastore
    //and the index key to id as well as the function to get data fron
    //the datastore (ensure column and index are loaded beofre carrying
    //out the function)
    setupLinks(dataStore,index,func){
        if (! this.browser){
            setTimeout(()=>this.setupLinks(dataStore,index,func),100)
        }
       
        this.dataLink= {
            dataStore:dataStore,
            index:index,
            getDataFunction:func
        }
        if (!this.bamscatrack){
            return;
        }
        const cat = this.config.cluster_reads;
      
        dataStore.addListener("gb_"+this.config.id,(type,data)=>{
            if (type==="filtered"){
                this.filterReads(data);
            }
            //has the data used to cluster reads changed
            //if so update the atac track
            else if (type==="data_changed"){
                const cl = this.config.cluster_reads;
                if (data.columns.indexOf(cl) !==-1){
                    this.changeClusters(cl,true)
                }
            }
        })
        func([cat],()=>{
            const ind = dataStore.getColumnIndex(index);
            const colors = dataStore.getColumnColors(cat);
            const col = dataStore.columnIndex[cat];
            this.bamscatrack.setCategories(cat,col.data,col.values,colors);
                //get basic dimension from the other datastore         
            this.cellDim= dataStore.getDimension("base_dimension");
          
            this.bamscatrack.addIndex(ind);
           
            if (this.config.color_by){
                this.colorByColumn(this.config.color_by);
            }
            
            //this.onDataHighlighted({indexes:[0]});
            /*const cf = dataStore.getColorFunction(cat);
            this.browser.setTrackColorFunction("_atac_bam_track",(feature)=>{
                const bc = ("test1#"+feature.tagBA.CB)+"";
                const ci = ind[bc];
                return cf(ci);
            })*/
            const g= this.config.genome_location;
            if (g){
                
                this.browser.update(g.chr,g.start,g.end,true);
            }
            else{
                this.onDataHighlighted({indexes:[0]});
            }  
        })
    }

    changeClusters(column,recalculateCats=false){
        this.config.cluster_reads=column;
        const ds = this.dataLink.dataStore;
        this.dataLink.getDataFunction([column],()=>{
            const colors = ds.getColumnColors(column);
            const col = ds.columnIndex[column];
            this.bamscatrack.setCategories(column,col.data,col.values,colors,true,recalculateCats);
            this.browser.update();
        });
           
    }
    themeChanged(){
        this.browser.repaint(true,true);
    }

    showHideDialog(){
        const b = this.browser;
        const controls = b.track_order.map(x=>{
            const c = b.tracks[x].config;
            return {
                type:"checkbox",
                id:c.track_id,
                label:c.short_label,
                current_value:c.hide?false:true
            }
        });
        new CustomDialog({
            title:"Show/Hide Tracks",
            controls:controls,
            doc:this.__doc__,
            buttons:[{
                text:"Update",
                method:vals=>{
                    for (let id in vals){
                        b.tracks[id].config.hide=!vals[id]
                    }
                    b.update();
                }
            }]
        })
      
    }

    getConfig(){
        const config= super.getConfig();
        const b = this.browser;
        config.tracks =   this.browser.getAllTrackConfigs();
        config.genome_location={
            chr:b.chr,
            start:b.start,
            end:b.end
        }
        return config;
    }

    filterReads(data){
        if (data===this.cellDim){
            return;
        }
        if (data==="all_removed"){
            this.bamscatrack.filterReads(null,null);
            this.resetButton.style.display="none";
            this.browser.removeHighlightedRegion("_filter");

        }
        else{
            this.bamscatrack.filterReads(this.dataLink.dataStore.filterArray,this.cellDim.filterArray);
        }
        this.browser.update();
    }
    //called if wig tracks associated with a column
    //obsolete???
    createColumnLinks(dataStore,columns,func){
        this.dataLink={
            dataStore:dataStore,
            columns:columns,
            getDataFunction:func
        }
        //only interested in the first column
        const col= columns[0].col;
        const vc= this.dataStore.getValueToColor(col);
        //color the wig tracks properly
        for (let v in vc){
            const tr =this.browser.tracks[`${col}|${v}`];
            if (tr){
                tr.config.color= vc[v];
            }
        }
       // this.onDataHighlighted({indexes:[0]});
    }

    //called when a feature is highlighted 
    onDataHighlighted(data){
        if (data.source===this){
            return;
        }       
        const p = this.config.param;
        const vm = this.config.view_margins;
        const o = this.dataStore.getRowAsObject(data.indexes[0],p);
        //some basic checks
        let st = o[p[2]]>o[p[1]]?o[p[1]]:o[p[2]];
        let en = o[p[2]]>o[p[1]]?o[p[2]]:o[p[1]];
        const rowData = data.data;

       
        if (rowData && rowData.phase_data){
            for (let id of this.browser.track_order){
                let pd = rowData.phase_data[id];           
                if (pd){
                    if (pd.ratio[0]==0 && pd.ratio[1]==0){
                        this.browser.setTrackAttribute(id,"phase_data",null)
                    } 
                    else{       
                        this.browser.setTrackAttribute(id,"phase_data",
                        {
                            ratio:pd.ratio[0]/(pd.ratio[0]+pd.ratio[1]),
                            st,
                            en
                        });
                    }
                }
                else{
                    this.browser.setTrackAttribute(id,"phase_data",null)
                }
            }
        }
        let margin =1000;
        if  (vm.type==="percentage"){
            en =  st===en?st+1:en;
            margin =  Math.round(en-st)*(vm.value/100);
        }
        else{
            margin = vm.value;
        }
        const fcp =this.config.feature_present_column;
        if (fcp){
            const ids =  this.dataStore.getRowText(data.indexes[0],fcp).split(", ");
            for (let id of ids){
                this.browser.setTrackAttribute(id,"color","#35de26")
            }
            const not =this.dataStore.getColumnValues(fcp).slice(0).filter(x=>ids.indexOf(x) ===-1);
            for (let id of not){
                this.browser.setTrackAttribute(id,"color","#939c92")
            }
            const new_order = ids.concat(not)
            const track_order=[];
            let index=0;
            for (let id of this.browser.track_order){
                if (new_order.indexOf(id) !==-1){
                    track_order.push(new_order[index]);
                    index++
                }
                else{
                    track_order.push(id)
                }
            }
            this.browser.track_order=track_order;
        }

        
        this.browser.update(o[p[0]],st-margin,en+margin) 
    }

    getImage(callback,type){
        if (type=="svg"){
            this.browser.repaint(true,false,callback)
        }
        else{
            callback(this.browser.canvas);
        }

    }

    configure_tracks(){
        this.dataSrore
    }


    getColorOptions(){		
		return {
			colorby:"all",
			has_default_color:true
		}    
	}

    _calculatePosition(text){
		text=text.replace(/,/g,"");
	
		let arr = text.split(":");
		let chr = null;
		let pos = null;
		if (arr.length===1){
			chr = this.browser.getPosition().chr;
			pos=arr[0]
		}
		else{
			chr =arr[0];
			pos=arr[1];
		}
		let arr2= pos.split("-");
		return ({chr:chr,start:parseInt(arr2[0]),end:parseInt(arr2[1])});
	}


    colorByColumn(column,update=true){
        const colorFunc =  this.getColorFunction(column);
        this.browser.setTrackColorFunction("_base_track",(feature)=>{
            return colorFunc(feature.data[0]);
        })
        if (update){
            this.browser.update();
        }
        
    }

    setLabelFunction(column){
        if (!column){
            this.browser.setTrackLabelFunction("_base_track",null);
            delete this.config.feature_label;
        }
        else{
            this.config.feature_label=column;
            this.browser.setTrackLabelFunction("_base_track",
            feature=>this.dataStore.getRowText(feature.data[0],column)     
            )
        }
    }

    onDataFiltered(data){
        if (!this.browser){
            setTimeout(()=>this.onDataFiltered(data),100);
            return;
        }
        if (!this.dataStore.isFiltered()){
            this.browser.setTrackFeatureFilter("_base_track",null);      
        }
        else{
            this.browser.setTrackFeatureFilter("_base_track",(feature)=>{
                if (this.dataStore.filterArray[feature.data[0]] ===0 ){
                    return true;
                }
                return false;
            });
        }
        
        this.browser.update();
    }

    removeFilter(){
        this.cellDim.removeFilter();
        this.browser.removeHighlightedRegion("_filter");
        this.browser.update();
    }

    changeBaseDocument(doc){
        this.browser.closeAllDialogs();
        super.changeBaseDocument(doc);
        this.browser.__doc__=doc;   
    }


    setSize(x,y){
        super.setSize(x,y);
        this.browser.setSize();
        this.browser.update();
    }

    remove(notify=true){
        if (this.cellDim){
            this.cellDim.destroy(notify);
            this.dataLink.dataStore.removeListener("gb_"+this.config.id);
        }
        super.remove();
    }

    getSettings(){
        const settings = super.getSettings();
        const cols = this.dataStore.getColumnList();
        const c = this.config;
        cols.push({name:"None",field:"_none"});
        if (this.bamscatrack){
            const cats= this.dataLink.dataStore.getColumnList("text");
            settings.push({
                label:"Group By",
                type:"dropdown",
                values:[cats,"name","field"],
                current_value:c.cluster_reads,
                    func:(x)=>{
                        this.changeClusters(x);
                    }
            });
        }
        settings.push({
            label:"Feature Label",
            type:"dropdown",
            values:[cols,"name","field"],
            current_value:c.feature_label || "_none",
                func:(x)=>{
                    x = x==="_none"?null:x;
                    this.setLabelFunction(x);
                    this.browser.update()
                }
        });
        settings.push({
            label:"View margin type",
            type:"radiobuttons",
            choices:[
                ["% of feature length","percentage"],
                ["absolute (bp)","absolute"]
            ],
            current_value:c.view_margins.type,
            func:(x)=>{
                c.view_margins.type=x;
                const d = this.dataStore.getHighlightedData();
                if (d){
                    this.onDataHighlighted({indexes:[d]});
                }
            }
        });
        settings.push({
            label:"View margin length",
            type:"text",
            current_value:c.view_margins.value,
            only_update_on_enter:true,
            func:(x)=>{
                x= parseInt(x);
                x= isNaN(x)?c.view_margins.type==="percentage"?20:1000:x;
                c.view_margins.value=x;
                const d = this.dataStore.getHighlightedData();
                if (d){
                    this.onDataHighlighted({indexes:[d]});
                }
            }
        });
        return settings;
    }
}

BaseChart.types["genome_browser"]={
    "class":GenomeBrowser,
    name:"Genome Browser",
    methodsUsingColumns:["setLabelFunction"],
    configEntriesUsingColumns:["feature_label","color_wig_tracks_by","feature_present_column"],

    params:[],
    required:["genome_browser"],
    init:(config,dataSource)=>{
        //set the chr,start,finish columns
        const gb = dataSource.genome_browser;
        config.param= gb.location_fields.slice(0);
        //legacy
        if (gb.feature_label){
            config.feature_label=gb.feature_label;
        }
        if (gb.default_parameters){
            Object.assign(config,gb.default_parameters);
        }
        
        //set the base track corresponding to the datastore
        const df= { 
            short_label:gb.default_track.label,
            url: gb.default_track.url,
            track_id:"_base_track",
            decode_function:"generic",
            displayMode:"EXPANDED"
        };
        config.tracks=[df];
        if (gb.default_track_parameters){
            Object.assign(df,gb.default_track_parameters)
        }
        //add the atac track if specified
        if (gb.atac_bam_track){
            //field to cluster reads
            config.cluster_reads=gb.atac_bam_track.cluster_reads;
            config.tracks.push({
                short_label:"Coverage",
                height:400,
                track_id:"_atac_bam_track",
                url:gb.atac_bam_track.url,
                type:"bam_sca_track"
            })
        }
        //add any default tracks
        config.default_track = gb.default_track.url;
        if (gb.default_tracks){
            for (let t of gb.default_tracks){
                config.tracks.push(JSON.parse(JSON.stringify(t)));
            }         
        }
    }
}

export default GenomeBrowser;