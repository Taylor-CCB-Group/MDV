
import BaseChart from "../BaseChart";
import { g } from  "../../lib/utils";
import SVCircosPlotComponent from "./SVCircosPlotComponent";

import { BaseReactChart } from "../../react/components/BaseReactChart";
import { action } from "mobx";

class SVCircosPlot extends BaseReactChart {
    constructor(dataStore,div,config) {
        if (!config.color_by){
            config.color_by=config.param[4];
        }
        config.flag =false;
        super(dataStore, div, config, SVCircosPlotComponent);
    }
    getColorOptions() {
        return {
            colorby: "all",
            has_default_color: false,
        };
    }
    colorByColumn(column){
        action(()=>{
            this.config.color_by= column;
        })();
      
    }
      getSettings(){
              const settings = super.getSettings();
              const c = this.config;
              return [
                ...settings,
                g({
                    //set whether filters are cleared on Reset All
                    type: "check",
                    current_value: c.flag || false,
                    label: "My Flag",
                    func: (x) => {
                       c.flag = x;
                    }
                })
              ]
        }

}


BaseChart.types["sv_circos_plot"] = {
    name: "SV Circos Plot",
    class: SVCircosPlot,
    //check the datasource can support svs
    required:(ds)=>{
        return ds.genome?.svs
    },
    //add the core sv columns to params
    init:(config,dataSource)=>{
        const cols= dataSource.genome?.svs?.sv_columns;
        if (!cols){
            throw Error("Trying to add SVCircosPlot, but No SV columns found in datasource genome metadata")
        }
        config.param=[cols["pos1"],cols["chr1"],cols["pos2"],cols["chr2"],cols["svtype"]];
    },
    //params are automatically set from genome metadata
    params: []
}

export default SVCircosPlot;