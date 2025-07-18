
import React from "react";
import { observer } from "mobx-react-lite";
import { registerChart } from "../chartTypes";
import schema from "./schema";
import {type MDVChartProps, useColumn} from "../MDVChart";



const HistogramChart = observer(({ chartInstance,datastore,filter }: MDVChartProps): React.ReactElement => {
    const x  = useColumn(chartInstance.x);
    function toggle(){
        if (!chartInstance.x.link) {
            chartInstance.x.setLink("genes","gs");
        }
        else{
            chartInstance.x.removeLink();
        }
        
    }
    return (
        <div>
            <h2>count:{datastore.filterSize}</h2>
            <h2>{x?x.name:"waiting"}</h2>
        
            <button 
                onClick={toggle}
            >
                {chartInstance.x.link ? `Linked to ${chartInstance.x.link.datastore} (${chartInstance.x.link.subgroup})` : "Set Link"}  
            </button>
      
        </div>
    );
});

registerChart("histogram", {ChartComponent:HistogramChart,schema});

export default HistogramChart;



