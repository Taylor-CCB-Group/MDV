import { DataModel } from "../table/DataModel.js";
import { Deck } from '@deck.gl/core/typed';
import BaseChart from "./BaseChart.js";
import { createEl } from "../utilities/Elements.js";
import { ImageArray } from "../webgl/ImageArray.js";
//import { ScatterplotLayer } from 'deck.gl/typed'; // -no, using ScatterplotExLayer
import { ScatterplotExLayer, ImageArrayDeckExtension } from '../webgl/ImageArrayDeckExtension.js';
import { OrbitView } from "deck.gl/typed";


class ImageScatterChart extends BaseChart {
    canvas: HTMLCanvasElement;
    imageArray: ImageArray;
    deck: Deck;
    dataModel: DataModel;
    constructor(dataStore, div, config) {
        super(dataStore, div, config);
        const canvas = this.canvas = createEl("canvas", {}, this.contentDiv);
        //const gl = canvas.getContext("webgl2"); // do we need to take care of disposing resources as well?
        const { base_url, image_key, texture_size } = config.images;
        this.dataModel = new DataModel(dataStore, { autoUpdate: false });
        this.imageArray = new ImageArray(dataStore, canvas, this.dataModel, {
            base_url,
            image_type: "png",
            image_key,
            width: texture_size,
            height: texture_size,
        });
        this.imageArray.onProgress = (n) => this.doDeck();
        const layers = this.getLayers();
        const view = new OrbitView();
        this.deck = new Deck({
            canvas,
            layers,
            views: [view],
            initialViewState: {
                // target: [0.5, 0.5, 0.5],
                // zoom: 4, //0 means "one pixel is one unit", 1 scales by 2
            },
            // glOptions: {},
            // parameters: {},
        });
    }
    doDeck() {
        const layers = this.getLayers();
        this.deck.setProps({ layers });
    }
    getLayers() {
        const { param } = this.config;
        const {columnIndex} = this.dataStore;
        type Column = {data: Float32Array, minMax: [number, number]}
        
        const cx = columnIndex[param[0]] as Column;
        const cy = columnIndex[param[1]] as Column;
        const cz = columnIndex[param[2]] as Column;
        function n(col, i) {
            const {minMax} = col;
            return 200*(col.data[i] - minMax[0]) / (minMax[1] - minMax[0]);
        }
        const data = new Array(cx.data.length).fill(0).map((_, i) => i);

        const {imageArray} = this;
        // const {getImageAspect, getImageIndex} = this.imageArray;// need to bind this
        type K = number;
        return [
            new ScatterplotExLayer({
                data,
                // baseRadius: 10000,
                billboard: true,
                getImageIndex: (i: K) => imageArray.getImageIndex(i),
                getImageAspect: (i: K) => imageArray.getImageAspect(i),
                // getPosition: (i: K) => [cx[i], cy[i], cz[i]],
                getPosition: (i: K) => [n(cx, i), n(cy, i), n(cz, i)],
                getRadius: (i: K) => 20,
                getFillColor: (i: K) => [255, 255, 255],
                imageArray, //XXX: deck doesn't like mutating props like this
                extensions: [new ImageArrayDeckExtension()]
            })
        ];
    }

    onDataFiltered() {
        this.dataModel.updateModel();
        this.doDeck();
    }
}

BaseChart.types["ImageScatterChart"] = {
    class: ImageScatterChart,
    name: "Image Scatter Plot",
    required: ["images"],
    methodsUsingColumns: ["onDataFiltered", "getLayers"],
    configEntriesUsingColumns: ["image_key"],

    init: (config, dataSource, extraControls) => {
        //get the available images
        const i = dataSource.images[extraControls.image_set];
        console.log('ImageScatterChart param', config.param);
        //set the base url and type
        config.images = {
            base_url: i.base_url,
            type: "png", //todo: allow this to be specified
            image_key: i.key_column, //nb ImageTableChart has this as config.param
            texture_size: extraControls.texture_size,
        };
        //allows configEntriesUsingColumns to work without this being a param.
        config.image_key = i.key_column;
    },
    extra_controls: (dataSource) => {
        const imageSets = [];
        for (let iname in dataSource.images) {
            imageSets.push({ name: iname, value: iname })
        }
        console.log('imageSets', imageSets);
        const sortableColumns = dataSource.getLoadedColumns().map(c => ({ name: c, value: c }));
        const imageSizes = [32, 64, 128, 256, 512, 1024].map(s => ({ name: s, value: s }));
        return [
            //drop down of available image sets
            {
                type: "dropdown",
                name: "image_set",
                label: "Image Set",
                values: imageSets
            },
            {
                type: "dropdown",
                name: "image_title",
                label: "Tooltip",
                values: sortableColumns
            },
            {
                type: "dropdown",
                name: "texture_size",
                label: "Texture Size",
                values: imageSizes
            },
        ];
    },
    params: [
        {
            type: "number",
            name: "X axis"
        },
        {
            type: "number",
            name: "Y axis"
        },
        {
            type: "number",
            name: "Z axis"
        },
        // ... some params should be optional
        // {
        //     type: "number",
        //     name: "radius"
        // },
    ]
};

export default ImageScatterChart;