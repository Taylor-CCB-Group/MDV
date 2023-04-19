import "microtip/microtip.css";
import "../src/css/fontawesome-5.15.3/all.css";
import 'nouislider/dist/nouislider.min.css'
import "../src/utilities/css/ContextMenu.css";
import "../src/charts/css/charts.css";
import "../src/webgl/css/wgl2di.css";
import "../src/table/css/slickgrid.css";

import ChartManager from "../src/charts/ChartManager";

import { Matrix4, Vector3 } from 'math.gl';
import { Plane } from '@math.gl/culling';
import { BaseDialog } from "../src/utilities/Dialog";
import { createEl } from "../src/utilities/Elements";
import { DataModel } from "../src/table/DataModel";
import { DataColumn } from "../src/charts/charts";


// let's add some mock data that looks a bit more like what we expect to see...
const size = 100000, numSlices = 5;
// these correspond to the size of the T1-Head sample.
const sizeX = 256, sizeY = 256, sizeZ = 129;

const goCellularTextVals = [
    "apicolateral plasma membrane",
    "nucleoplasm",
    "histone acetyltransferase complex",
    "trans-Golgi network transport vesicle",
    "nuclear outer membrane",
    "fusome",
    "presynaptic periactive zone",
    "nucleolus",
    "type Ib terminal bouton",
    "nBAF complex",
    "chaperonin-containing T-complex",
    "integrin complex",
    "Wnt signalosome",
    "female germline ring canal inner rim",
    "endoplasmic reticulum membrane",
    "spanning component of plasma membrane",
    "extrinsic component of mitochondrial outer membrane",
    "nuclear envelope",
    "axon",
    "actin cytoskeleton",
    "ruffle",
    "protein-containing complex",
    "apical plasma membrane",
    "meiotic nuclear membrane microtubule tethering complex",
    "nuclear periphery",
    "astral microtubule",
    "focal adhesion",
    "precatalytic spliceosome",
    "proton-transporting ATP synthase complex",
    "cytoplasmic vesicle",
    "Mpp10 complex",
    "U2 snRNP",
    "cytoskeleton of presynaptic active zone",
    "cytoplasmic stress granule",
    "lysosome",
    "apical part of cell",
    "neuronal cell body",
    "lysosomal membrane",
    "SWI\/SNF complex",
    "endoplasmic reticulum exit site",
    "receptor complex",
    "growth cone",
    "pleated septate junction",
    "extrinsic component of plasma membrane",
    "sodium:potassium-exchanging ATPase complex",
    "heterochromatin",
    "dystrophin-associated glycoprotein complex",
    "germline ring canal",
    "Rb-E2F complex",
    "perivitelline space",
    "tricellular tight junction",
    "mitotic spindle pole",
    "peroxisome",
    "vesicle membrane",
    "somatic ring canal",
    "extracellular region",
    "endocytic vesicle",
    "investment cone",
    "presynaptic active zone",
    "U1 snRNP",
    "membrane",
    "filopodium",
    "catalytic core F(1)",
    "Wnt-Frizzled-LRP5\/6 complex",
    "catalytic core",
    "A band",
    "ASAP complex",
    "sarcomere",
    "beta-catenin-TCF complex",
    "perinuclear region of cytoplasm",
    "axonal growth cone",
    "postsynaptic density membrane",
    "cytoneme"
]

const mockDataSources = [{
    size,
    name: "test", //whatever...
    columns: [
        {
            name: "cell_id",
            datatype: "unique"
        },
        {
            name: "slice_id",
            datatype: "text"
        },
        {
            name: "mockX",
            datatype: "double"
        },
        {
            name: "mockY",
            datatype: "double"
        },
        {
            name: "mockZ",
            datatype: "double"
        },
        {
            name: "classifier",
            datatype: "text"
        },
        {
            name: "GO_cellular_component",
            datatype: "multitext",
            values: goCellularTextVals
        }
    ],
    columnGroups: [{
        name: "spatial data",
        columns: ["mockX", "mockY", "mockZ"]
    }]
}];
const U = undefined as any; //YOLO
const mockSlices = new Array(numSlices).fill(U).map((_, slice_id) => {
    // make a bunch of points, all on a random plane...
    const r = (v=1) => v*Math.random();
    // const normal = new Vector3([r(2)-1, r(2)-1, r(2)-1]);
    const normal = new Vector3([r(0.5)-0.25, -1, r(0.5)-0.25]);
    normal.normalize();
    const d = slice_id * sizeY / numSlices; //r(32);
    const plane = new Plane(normal, d*0.8);
    const p = (_, i) => {
        const p2d = [r(256), 0, r(256)];
        const p = plane.projectPointOntoPlane(p2d);
        const mockX = p[0], mockY = p[1], mockZ = p[2];
        const classifier = "Type " + Math.round(r(3));
        const w = Math.pow(1-r(), 3) * (goCellularTextVals.length-1);
        const GO_cellular_component = goCellularTextVals[Math.floor(w)];
        /// not normalized wrt slice_id / plane... I suppose I need to figure out how dataSources relate
        const cell_id = slice_id + '-' + i;
        return {cell_id, slice_id, mockX, mockY, mockZ, classifier, plane, GO_cellular_component};
    }
    const n = size / numSlices;
    const points = new Array(n).fill(U).map(p);
    return {
        slice_id,
        plane,
        points
    }
});

// flatten mockSlices...
const mockDataRows = [];
for (const s of mockSlices) {
    mockDataRows.push(...s.points);
}
//const mockDataCols = new Map();

const mockDataLoaderFn = async (columns, dataSource, size) => {
    const dataList: {field: string, data}[] = [];
    for (const column of columns) {
        const data = mockDataRows.map(r => r[column.name]);
        dataList.push({ field: column.name, data });
    }
    return dataList;
}

const mockDataLoader = {
    function: mockDataLoaderFn
}

// could we add a 'select slice' option, that would 
// - select the relevant points
// - add an image-layer to the 3d view with an image & orientation corresponding to the slice
// - set a clipping plane on the 3d view

class SliceDialog extends BaseDialog {
    constructor(chartManager: ChartManager) {
        super({
            title: "Select slice",
            width: 300,
            height: 200,
            columns: 1,
        }, {});
        const dataStore = chartManager.getDataSource('test').dataStore;
        const sliceColumn = dataStore.getColumn('slice_id') as DataColumn<'text'>;
        // create a dropdown with the slice ids
        const sliceDropdown = createEl('select', {}, this.columns[0]);
        for (const s of mockSlices) {
            const option = document.createElement('option');
            sliceDropdown.appendChild(option);
            option.value = '' + s.slice_id;
            option.innerText = 'slice ' + s.slice_id;
        }
        sliceDropdown.addEventListener('change', (e) => {
            const slice_id = Number.parseInt(sliceDropdown.value);
            const slice = mockSlices.find(s => s.slice_id === slice_id);
            //chartManager.getDataSource('test').selectRows(r => r.slice_id === slice_id);
            // how do we filter the data?
                            
        });
    }
    // init() {
    // }
    // open(callback) {
    // }
}

 //config - single view
const config={
    gridstack: true,
    only_view:{
        //the initial charts to load
        initialCharts:{
            "test":[
            {
                type: "viv_volume_scatter_view",
                title: "T1-head volume view",
                param: ["mockX", "mockY", "mockZ"],
                options: {url: "/data/t1-head-imj.ome.tiff"},
                size:[400,600],
                position:[10,10]
            },
            {//TODO figure out why 'field colors' isn't working
                type: "table_chart",
                title: "table",
                id: "table",
                param: ["cell_id", "slice_id", "classifier"],
                size: [400, 600],
                position: [440, 10]
            },
            {
                type: "wgl_scatter_plot",
                title: "mockX x mockY",
                id: "scatter",
                param: ["mockX", "mockY"],
                size: [400, 300],
                position: [10, 620],
                color_by: "slice_id"
            }
            // {
            //     type:"row_chart",
            //     param:"Clusters",
            //     id:"my-clusters",
            //     size:[200,390],
            //     position:[420,220],
            // }
        
            ]
        }
    }
};
//create the app with the above parameters
// const cm = new ChartManager("app1", dataSources, dataLoader, config);
const cm = new ChartManager("app1", mockDataSources, mockDataLoader, config);
export default cm;