import "microtip/microtip.css";
import "../css/fontawesome-5.15.3/all.css";
import 'nouislider/dist/nouislider.min.css'
import "../utilities/css/ContextMenu.css";
import "../charts/css/charts.css";
import "../webgl/css/wgl2di.css";
import "../table/css/slickgrid.css";

import {getRandomDataStore} from "../dev/MockDataStore.js";
import { createEl } from "../utilities/Elements";
import  RingChart  from "../charts/RingChart";
import  RowChart  from "../charts/RowChart";
import  WGLScatterPlot  from "../charts/WGLScatterPlot";
import  WGL3DScatterPlot  from "../charts/WGL3DScatterPlot";
import TableChart from "../charts/TableChart";
import  HistogramChart  from "../charts/HistogramChart";


const ds2= getRandomDataStore(10000000);
const ap= document.getElementById("app1");

const g1 = createEl("div",{classes:["holder"]},ap);
new RingChart(ds2,g1,{param:"colors"});
const g2 = createEl("div",{classes:["holder"]},ap);
new RowChart(ds2,g2,{param:["colors"]});
const g3 = createEl("div",{classes:["holder"]},ap);
new WGLScatterPlot(ds2,g3,{param:["x1","x2"],color_by:"x1",color_legend:{display:true,pos:[50,50]}});
const g4 = createEl("div",{classes:["holder"]},ap);
new WGL3DScatterPlot(ds2,g4,{param:["x1","x2","x3"],color_by:"colors"});
const g5 = createEl("div",{classes:["holder"]},ap);
new TableChart(ds2,g5,{param:["x1","x2","x3","colors"]});
const g6 = createEl("div",{classes:["holder"]},ap);
new HistogramChart(ds2,g6,{param:["x3"]});












