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


// let's add some mock data that looks a bit more like what we expect to see...
const size = 10000, numSlices = 5;
// these correspond to the size of the T1-Head sample.
const sizeX = 256, sizeY = 256, sizeZ = 129;

const mockDataSources = [{
    size,
    name: "cells", //whatever...
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
        }
    ],
    columnGroups: [{
        name: "spatial data",
        columns: ["mockX", "mockY", "mockZ"]
    }]
}];

const mockSlices = new Array(numSlices).fill().map((_, slice_id) => {
    // make a bunch of points, all on a random plane...
    const r = (v=1) => v*Math.random();
    const normal = new Vector3([r(2)-1, r(2)-1, r(2)-1]);
    normal.normalize();
    const d = r(256);
    const plane = new Plane(normal, d);
    const p = (_, i) => {
        const p2d = [r(256), r(256), 0];
        const p = plane.projectPointOntoPlane(p2d);
        const mockX = p[0], mockY = p[1], mockZ = p[2];
        const classifier = "Type " + Math.round(r(3));
        /// not normalized wrt slice_id / plane... I suppose I need to figure out how dataSources relate
        const cell_id = slice_id + '-' + i;
        return {cell_id, slice_id, mockX, mockY, mockZ, classifier, plane};
    }
    const n = 10000 / 5;
    const points = new Array(n).fill().map(p);
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
    const dataList = [];
    for (const column of columns) {
        const data = mockDataRows.map(r => r[column.name]);
        dataList.push({ field: column.name, data });
    }
    return dataList;
}

const mockDataLoader = {
    function: mockDataLoaderFn
}


 //config - single view
const config={
    only_view:{
        //the initial charts to load
        initialCharts:{
            "cells":[
            {
                type: "viv_volume_view",
                title: "T1-head volume view",
                param: [],
                options: {url: "/data/t1-head-imj.ome.tiff"},
                size:[400,600],
                position:[10,10]
            },
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