import BaseChart from "../BaseChart";
import ReactDOM from "react-dom/client";
import { createElement } from "react";
import {test_config} from "./config_utils";
import { createEl } from "../../utilities/Elements";
//import {createViewState,JBrowseApp} from '@jbrowse/react-app2'
//import makeWorkerInstance from '@jbrowse/react-app2/esm/makeWorkerInstance'

class  JBrowse extends BaseChart {
    constructor(dataStore,div,config) {
        // Extract jbrowse_config and prevent MobX double-wrapping
        let jbrowse_config = config.jbrowse_config;
        config.jbrowse_config = undefined;
        super(dataStore, div, config);

        // Create the browser container div
        this.browserDiv = createEl("div", {
            styles: {
                height: "100%",
                position: "relative",
                overflowY: "auto"
            }
        }, this.contentDiv);

        // Use provided jbrowse_config or create a default one
        this.jbrowse_config = jbrowse_config || this.createConfig();

        // Initialize JBrowse asynchronously
        this.init();
    }

    createConfig(){
        // ToDO: customize config based on dataStore or other parameters
       return test_config;
    }

    /**
     * Asynchronously loads JBrowse modules and renders the genome browser.
     */
    async init() {
        try {
            // Dynamically import JBrowse modules
            const { createViewState, JBrowseApp } = await import("@jbrowse/react-app2");
            const { default: makeWorkerInstance } = await import("@jbrowse/react-app2/esm/makeWorkerInstance");

            // Create the JBrowse view state
            this.viewState = createViewState({
                config: this.jbrowse_config,
                makeWorkerInstance,
            });
            //need to change  if breakpair view
            this.mainView= this.viewState.session.views[0];

            // Render the JBrowse React app into the browserDiv
            this.root = ReactDOM.createRoot(this.browserDiv);
            this.root.render(
                createElement(JBrowseApp, { viewState: this.viewState })
            );
        } catch (error) {
            // Handle dynamic import or render errors gracefully
            console.error("Failed to initialize JBrowse:", error);
            this.browserDiv.innerHTML = "<div style='color:red'>Failed to load genome browser.</div>";
        }
    }

    onDataHighlighted(data){
        console.log(data);
        const index = data.indexes[0];
        const p = this.config.param;
        if (this.dataStore.genome?.svs){
            //get all the sv information
            const chr1 = this.dataStore.getRowText(index,p[0]);
            const pos1 = this.dataStore.getRowText(index,p[1]);
            const chr2 = this.dataStore.getRowText(index,p[2]);
            const pos2 = this.dataStore.getRowText(index,p[3]);
            //on different chromosomes just show the first position
            //breakpoint view can handle the rest
            if (chr1 !== chr2 || pos2-pos1 > 20000){
                this.setLocation(chr1,pos1-500,pos1+500);
            }
            //if on same chromosome and not too far apart show the whole region
            else{
                this.setLocation(chr1,pos1-500,pos2+500);
            }
        }
    }

    setLocation(chr,start,end){
        if (this.mainView){
            this.mainView.navToLocString(`${chr}:${start}-${end}`);
        }
    }
  
    remove(notify) {
        if (this.root) {
            this.root.unmount();
        }
        super.remove(notify);
    }

    getConfig() {
        const c = super.getConfig();
        let jbrowse_config = undefined;

        // If viewState and session exist, serialize the session
        if (this.viewState && this.viewState.session) {
            jbrowse_config = {
                ...this.jbrowse_config,
                defaultSession: this.viewState.session.toJSON()
            };
            c.jbrowse_config = jbrowse_config;
        } else {
            // Otherwise, just use the initial config
            c.jbrowse_config = this.jbrowse_config;
        }
        return c;
    }
}

BaseChart.types["jbrowse_genome_viewer"] = {
    name: "Genome Browser (JBrowse)",
    class: JBrowse,
    //check the datasource can support svs
    required:["genome"],
    //add the core sv columns to params
    init:(config,dataSource)=>{
        const svs = dataSource.genome?.svs;
        //if this genome datasource represents svs
        if (svs){
             const cols= svs.sv_columns;
            if (!cols){
                throw Error("Trying to add SVCircosPlot, but No SV columns found in datasource genome metadata")
            }
        config.param=[cols["pos1"],cols["chr1"],cols["pos2"],cols["chr2"]];
        }
       
    },
    //params are automatically set from genome metadata
    params: []
}

export default JBrowse