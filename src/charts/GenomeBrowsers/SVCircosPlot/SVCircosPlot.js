
import BaseChart from "../../BaseChart";
import SVCircosPlotComponent from "./SVCircosPlotComponent";
import { BaseReactChart } from "../../../react/components/BaseReactChart";




class SVCircosPlot extends BaseReactChart {
    colorVersion = 0;
    constructor(dataStore,div,config) {
        //check the datasource has svs and chromosomes
        //could this be done with zod schema ? 
        if (!dataStore.genome || dataStore.genome.type !== "sv" || !dataStore.genome.columns) {
            throw new Error("DataStore genome must have svs and chromosomes for SVCircosPlot");
        } 
        //initially color by sv type
        if (!config.color_by){
            config.color_by=dataStore.genome.columns.svtype;
        }
        //all this could be automatic bases on zod scheama
        config.trim_color_scale = config.trim_color_scale  ??  "none";
        config.fallbackOnZero = config.fallbackOnZero ?? false;
        config.hideMissing = config.hideMissing ?? false;
        config.log_color_scale = config.log_color_scale ?? false;
        super(dataStore, div, config, SVCircosPlotComponent);
    }
    getColorOptions() {
        return {
            colorby: "all",
            has_default_color: false,
        };
    }
    colorByColumn(column){
        this.config.color_by = column;
    }
    
}


BaseChart.types["sv_circos_plot"] = {
    name: "SV Circos Plot",
    class: SVCircosPlot,
    //check the datasource can support svs
    required:(ds)=>{
        return ds.genome && ds.genome.type==="sv" && ds.genome.chromosomes && ds.genome.columns
    },
}

export default SVCircosPlot;
