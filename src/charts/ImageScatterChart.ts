import { DataModel } from "../table/DataModel.js";
import { Deck } from '@deck.gl/core/typed';
import BaseChart from "./BaseChart.js";
import { createEl } from "../utilities/Elements.js";
import { ImageArray } from "../webgl/ImageArray.js";
//import { ScatterplotLayer } from 'deck.gl/typed'; // -no, using ScatterplotExLayer
import { ScatterplotExLayer, ImageArrayDeckExtension } from '../webgl/ImageArrayDeckExtension.js';
import { OrbitView } from "deck.gl/typed";

// not a definitive type, but marginally better than 'any', locally for now...
type Column = { data: Float32Array, minMax: [number, number] }
let nextID = 0;
class ImageScatterChart extends BaseChart {
    canvas: HTMLCanvasElement;
    imageArray: ImageArray;
    deck: Deck;
    dataModel: DataModel;
    progress = 0;
    billboard = true;
    size = 20;
    _layer: ScatterplotExLayer;
    id: number;
    constructor(dataStore, div, config) {
        super(dataStore, div, config);
        this.id = nextID++;
        const canvas = this.canvas = createEl("canvas", {}, this.contentDiv);
        //const gl = canvas.getContext("webgl2"); // do we need to take care of disposing resources as well?
        const { base_url, image_key, texture_size } = config.images;
        this.dataModel = new DataModel(dataStore, { autoUpdate: false });
        
        //---- PJT XXX::: this always trips me up... wasting too much time here...
        this.dataModel.setColumns([...config.param, image_key, config.image_title]);
        this.dataModel.updateModel();
        //----------------

        this.imageArray = new ImageArray(dataStore, canvas, this.dataModel, {
            base_url,
            image_type: "png",
            image_key,
            width: texture_size,
            height: texture_size,
        });
        this.imageArray.onProgress = (n) => {
            this.progress = n;
            this._layer.setNeedsUpdate();// attempt to force redraw
            this.doDeck();
        }
        const layers = this.getLayers();
        const view = new OrbitView();
        this.deck = new Deck({
            canvas,
            layers,
            views: [view],
            controller: true,
            initialViewState: {
                // if these are not set, there is an error when first using mouse-wheel to zoom
                target: [0, 0, 0], 
                zoom: 0, //0 means "one pixel is one unit", 1 scales by 2
            },
            getTooltip: (info) => {
                const {index, picked} = info;
                const titleColumn = this.config.image_title;
                // this.dataModel.updateModel();
                const text = this.dataModel.getItemField(index, titleColumn);
                return picked && {html: `<div>${titleColumn}: '${text}'</div>`,}
            },
            // glOptions: {},
            // parameters: {},
        });
    }
    doDeck() {
        // not necessary ATM... layers should update with updateTriggers
        // that's not working, but this also isn't helping...
        const layers = this.getLayers();
        this.deck.setProps({ layers });
    }
    getLayers() {
        const { param } = this.config;
        const {columnIndex} = this.dataStore;
        
        const cx = columnIndex[param[0]] as Column;
        const cy = columnIndex[param[1]] as Column;
        const cz = columnIndex[param[2]] as Column;
        function n(col: Column, i: number) {
            const {minMax} = col;
            return 200*(col.data[i] - minMax[0]) / (minMax[1] - minMax[0]) - 100;
        }
        
        const {length} = cx.data;
        const {progress, size, billboard} = this;
        
        // only length is needed - trying to force more updates, 
        // but for performance when working properly we should target updates more precisely
        const data = {length, progress, size, billboard}; //new Uint32Array(cx.data.length+1).fill(0).map((_, i) => i);
        //data[cx.data.length] = this.progress;//experiment wirh shallow difference to try to force update

        const {imageArray} = this;
        // const {getImageAspect, getImageIndex} = this.imageArray;// need to bind this
        type K = never; //data is not iterable, so ~never will be passed to the callbacks
        const layer = new ScatterplotExLayer({
            id: this.id,
            data,
            // radiusUnits: 'pixels', //default 'meters', also lineWidthUnits...
            // stroked: true, //TODO: make sure we can render properly with this
            billboard: this.billboard,
            pickable: true,
            getImageIndex: (i: K, {index}) => imageArray.getImageIndex(index),
            getImageAspect: (i: K, {index}) => imageArray.getImageAspect(index),
            getPosition: (_: K, {index, data, target}) => {
                const i = index as number;
                //[n(cx, i), n(cy, i), n(cz, i)] // say no to garbage
                target[0] = n(cx, i);
                target[1] = n(cy, i);
                target[2] = n(cz, i);
                return target;
            },
            getRadius: 1,
            radiusScale: this.size,
            getFillColor: [255, 255, 255],
            imageArray,
            updateTriggers: {
                //TODO: figure out how to get this to work...
                getImageAspect: [this.progress],
                // imageAspect: [this.progress],
            },
            extensions: [new ImageArrayDeckExtension()]
        });
        this._layer = layer;
        return [layer];
    }

    onDataFiltered() {
        this.dataModel.updateModel();
        this.doDeck();
    }

    getSettings() {
        return [...super.getSettings(),
            {
                name: "billboard",
                label: "Billboard",
                type: "check",
                current_value: this.billboard,
                func: (v) => {
                    this.billboard = v;
                    this.doDeck();
                }
            },
            {
                type: "slider",
                name: "size",
                label: "Size",
                current_value: this.size,
                min: 5,
                max: 100,
                step: 1,
                continuous: true,
                func: (v) => {
                    this.size = v;
                    this.doDeck();
                }
            },
        ]
    }
}

BaseChart.types["ImageScatterChart"] = {
    class: ImageScatterChart,
    name: "Image Scatter Plot",
    required: ["images"],
    methodsUsingColumns: ["onDataFiltered", "getLayers"],
    configEntriesUsingColumns: ["image_key", "image_title"],

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
                values: imageSizes,
                defaultVal: 256
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