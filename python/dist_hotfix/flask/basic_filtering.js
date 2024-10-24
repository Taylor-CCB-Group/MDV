


import "microtip/microtip.css";
import "../src/css/fontawesome-5.15.3/all.css";
import 'nouislider/dist/nouislider.min.css'
import "../src/utilities/css/ContextMenu.css";
import "../src/charts/css/charts.css";
import "../src/webgl/css/wgl2di.css";
import "../src/table/css/slickgrid.css";

import {getRandomDataStore} from "../src/dev/MockDataStore.js";
import { createEl } from "../src/utilities/Elements";
import  WGLScatterPlot  from "../src/charts/WGLScatterPlot";

const ap= document.getElementById("app1");
const div = createEl("div",{classes:["holder"]},ap);
const ds= getRandomDataStore(1000000);
new WGLScatterPlot(ds,div,{param:["x1","x2"],color_by:"x1",color_legend:{display:true,pos:[50,50]}});

//filter columns x1 and x2 both to >-25  and <25
//only center points visible in scatter plot
const cols  = ["x1","x2"];
const ranges={range1:[-25,25],range2:[-25,25]};
const r_dim = ds.getDimension("range_dimension");
r_dim.filter("filterSquare",cols,ranges)

