import "microtip/microtip.css";
import "../src/css/fontawesome-5.15.3/all.css";
import 'nouislider/dist/nouislider.min.css'
import "../src/utilities/css/ContextMenu.css";
import "../src/charts/css/charts.css";
import "../src/webgl/css/wgl2di.css";
import "../src/table/css/slickgrid.css";

import {getRandomDataStore, mockColumns, mockDataLoader, mockDataSource} from "../src/dev/MockDataStore.js";
import { createEl } from "../src/utilities/Elements";
import  RingChart  from "../src/charts/RingChart";
import  RowChart  from "../src/charts/RowChart";
import  WGLScatterPlot  from "../src/charts/WGLScatterPlot";
import  WGL3DScatterPlot  from "../src/charts/WGL3DScatterPlot";
import TableChart from "../src/charts/TableChart";
import  HistogramChart  from "../src/charts/HistogramChart";
import ChartManager from "../src/charts/ChartManager";

////////////////////
////////////////////
////////
//////// Original 10M_points example didn't use ChartManager...
//////// adapting so that it has a more full UI, using a dataLoader

//const ds2= getRandomDataStore(10000001);
const dataSource = mockDataSource("ds0", 10000001);
// const ap= document.getElementById("app1");

// const g1 = createEl("div",{classes:["holder"]},ap);
// new RingChart(ds2,g1,{param:"colors"});
// const g2 = createEl("div",{classes:["holder"]},ap);
// new RowChart(ds2,g2,{param:["colors"]});
// const g3 = createEl("div",{classes:["holder"]},ap);
// new WGLScatterPlot(ds2,g3,{param:["x1","x2"],color_by:"x1",color_legend:{display:true,pos:[50,50]}});
// const g4 = createEl("div",{classes:["holder"]},ap);
// new WGL3DScatterPlot(ds2,g4,{param:["x1","x2","x3"],color_by:"colors"});
// const g5 = createEl("div",{classes:["holder"]},ap);
// new TableChart(ds2,g5,{param:["x1","x2","x3","colors"]});
// const g6 = createEl("div",{classes:["holder"]},ap);
// new HistogramChart(ds2,g6,{param:["x3"]});


const config = {
    gridstack: true,
    only_view: {
        initialCharts: {"ds0": [
            { type: "ring_chart", param: ["colors"] },
            { type: "row_chart", param: ["colors"] },
            { type: "wgl_scatter_plot", param: ["x1", "x2"], color_by: "x1", color_legend: { display: true, pos: [50, 50] } },
            { type: "wgl_3d_scatter_plot", param: ["x1", "x2", "x3", "colors"], color_by: "colors" },
            { type: "wgl_3d_scatter_plot", param: ["x1", "x2", "x3", "colors"] },
            { type: "bar_chart", param: ["x3"]}
        ]}
    }
}

const mockLoader = mockDataLoader();

new ChartManager("app1", [dataSource], mockLoader, config);