import { DataModel } from "../table/DataModel.js";
import { Deck } from '@deck.gl/core/typed';
import {
    loadOmeTiff,
    DetailView,
    DETAIL_VIEW_ID,
    getChannelStats
} from '@hms-dbmi/viv';

import BaseChart from "./BaseChart.js";
import { createEl } from "../utilities/Elements.js";
import { ImageArray } from "../webgl/ImageArray";
import { OrthographicView, ScatterplotLayer } from "deck.gl/typed";
import Dimension from "../datastore/Dimension.js";
import { getProjectURL } from "../dataloaders/DataLoaderUtil.js";

// not a definitive type, but marginally better than 'any', locally for now...
type Column = { data: Float32Array, minMax: [number, number] }

type OMEConfig = {
    type: "ImageScatterTest"
    legend: string,
    param: string[],
    title: string,
    imageURL: string,
    color_by?: string,
    color_legend?: any,
}

let nextID = 0;
class VivScatterPlotNew extends BaseChart {
    canvas: HTMLCanvasElement;
    imageArray: ImageArray;
    deck: Deck;
    detailView: DetailView;
    tiffLoader: Awaited<ReturnType<typeof loadOmeTiff>>;
    dataModel: DataModel;
    progress = 0;
    billboard = true;
    size = 13;
    opacity = 255;
    saturation = 1;
    colorBy?: (index: number) => number[];
    id: number;
    constructor(dataStore, div, config: OMEConfig) {
        super(dataStore, div, config);
        console.warn('VivScatterPlotNew is only for testing, and will almost certainly be removed soon.');
        this.id = nextID++;
        const canvas = this.canvas = createEl("canvas", {}, this.contentDiv);
        this.dataModel = new DataModel(dataStore, { autoupdate: true });
        this.dataModel.setColumns(config.param);
        this.dataModel.updateModel();
        // in future, something slightly more sophisticated than this...
        loadOmeTiff(config.imageURL).then((tiff) => {
            this.tiffLoader = tiff;

            const width = tiff.metadata.Pixels.SizeX;
            const height = tiff.metadata.Pixels.SizeY;

            this.detailView = new DetailView({
                id: DETAIL_VIEW_ID + this.id,
                width,
                height,
            });

            // --- getDefaultChannelStats etc... refer to VivViewerMDV createLayers() &co ---
            
            const layers = this.updateDeck(); //...
            const view = new OrthographicView({});
            this.deck = new Deck({
                canvas,
                layers,
                views: [view],
                controller: true,
                initialViewState: {
                    // if these are not set, there is an error when first using mouse-wheel to zoom
                    target: [0, 0, 0],
                    zoom: -8, //0 means "one pixel is one unit", 1 scales by 2
                },
            });
            if (config.color_by) {
                this.colorByColumn(config.color_by);
            }
        });
    }
    updateDeck() {
        const { param } = this.config;
        const { columnIndex } = this.dataStore;

        const cx = columnIndex[param[0]] as Column;
        const cy = columnIndex[param[1]] as Column;
        
        const { data } = this.dataModel;
        type K = number;

        function n(col: Column, i: number) {
            const { minMax } = col;
            //TODO scaling options
            return (col.data[i] - minMax[0]) / (minMax[1] - minMax[0]);
        }
        // XXX: I'm not seeing anything as yet... spector.js isn't showing any xr-layer shader calls...
        const vivLayers = this.detailView.getLayers({
            viewStates: {
                id: `viv-${this.id}`,
                // target: [0, 0, 0],
                // zoom: -8,
            }, //TODO
            props: {
                loader: this.tiffLoader.data,
                // todo: un-hardcode this, sensible defaults & configurable options
                contrastLimits: [[98, 63307]],
                channelsVisible: [true],
                selections: [{z: 0, c: 3, t: 0}]
            }
        });

        const layer = new ScatterplotLayer({
            id: `scatter-${this.id}`,
            data,
            pickable: true,
            getPosition: (i: K, { target }) => {
                // say no to garbage (probably doesn't matter with generational GC & this being in 'nursery')
                //[n(cx, i), n(cy, i), n(cz, i)]
                target[0] = cx.data[i];
                target[1] = cy.data[i];
                target[2] = 0;//n(cz, i);
                return target as unknown as Float32Array; //🤮
            },
            getRadius: 1,
            radiusScale: this.size,
            // getFillColor: this.colorBy ? (i: K)=>[...this.colorBy(i), this.opacity] : [255, 255, 255, this.opacity],
            opacity: this.opacity ** 2,
            getFillColor: this.colorBy ? (i: K) => this.colorBy(i) : [255, 255, 255] as any,
            updateTriggers: {
                getFillColor: [this.colorBy, this.opacity],
            },
            // extensions: [new ImageArrayDeckExtension()]
        });
        const layers = [...vivLayers, layer];
        if (this.deck) this.deck.setProps({ layers });
        return layers;
    }

    onDataFiltered(dim: Dimension) {
        this.dataModel.updateModel();
        this.updateDeck();
    }

    colorByColumn(col: string) {
        this.colorBy = this.getColorFunction(col, true);
        this.updateDeck();
    }
    colorByDefault() {
        this.colorBy = null;
        this.updateDeck();
    }


    getColorOptions() {
        return {
            colorby: "all",
        }
    }

    getSettings() {
        return [...super.getSettings(),
        {
            type: "slider",
            name: "size",
            label: "Size",
            current_value: this.size,
            min: 1,
            max: 200,
            step: 1,
            continuous: true,
            func: (v) => {
                this.size = v;
                this.updateDeck();
            }
        },
        {
            type: "slider",
            name: "opacity",
            label: "Opacity",
            current_value: this.opacity,
            min: 0,
            max: 1,
            // step: 1,
            continuous: true,
            func: (v) => {
                this.opacity = v;
                this.updateDeck();
            }
        },
        ]
    }
}

// experiment with something similar to "viv_scatter_plot"
// how does this relate to "image_scatter_plot" / DensityScatterPlot?
// ---> actually, those are too complicated in various ways, and following the un-typed code through is a pain.
//      I'm going to make something with a much more basic option for imageURL... no use of 'regions' code for now.
BaseChart.types["VivScatterPlotNew"] = {
    class: VivScatterPlotNew,
    name: "Viv Scatter Plot (New)",
    // required: ds => ds.regions?.avivator,
    params: [
        {
            type: "number",
            name: "x",
        },
        {
            type: "number",
            name: "y",
        }
    ],
    extra_controls: (ds) => {
        return [
            {
                type: "string",
                name: "imageURL",
                label: "ome.tiff URL",
                defaultVal: ""
            }
        ]
    },
}

export default VivScatterPlotNew;