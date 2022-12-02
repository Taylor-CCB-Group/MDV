import BaseChart from "../charts/BaseChart.js";
import { ColorChannelDialog } from "../charts/VivScatterPlot.js";
import { createEl } from "../utilities/Elements.js";
import VivViewer from '../webgl/VivViewer.js';

class VivVolume extends BaseChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config);
        this.setupScatterplot();
        this.afterAppCreation();
        this.addMenuIcon("fas fa-palette", "Alter Channels").addEventListener("click", ()=>{
            if (this.dialog) {
                this.dialog.close();
                this.dialog = null;
            }
            else this.dialog = new ColorChannelDialog(this.viv);
        });
    }
    setSize(x, y) {
        super.setSize(x, y);
        if (this.viv) {
            const b = this._getContentDimensions();
            this.viv.setSize(b.width, b.height, {});
        }
    }
    //could just be in constructor in this case.
    afterAppCreation() {
        const {width, height} = this._getContentDimensions();
        this.vivCanvas = createEl("canvas", {
            width, height, 
            styles: {
                position: "absolute"
            }
        }, this.contentDiv);
        console.log('config in VivVolume', this.config);
        const vivConfig = {
            ///...
            url: this.config.options.url, //'/data/t1-head-imj.ome.tiff',
            use3d: true,
            scatterData: this.scatterData,
            getScatterFillColor: () => [100, 0, 0]
        }
        // not really what we want... basically ignored with use3d.
        const iv = {
            x_scale: 1, y_scale: 1, offset: 0
        }
        // we may not necessarily want to re-use VivViewer?
        // although some common features with channels etc.
        this.viv = new VivViewer(this.vivCanvas, vivConfig, iv);
    }

    setupScatterplot() {
        const {dataStore, config} = this;
        if (config.param.length === 0) return;
        //we'll need all the methods for handling data filtering etc...
        const ix = config.param[0];
        const iy = config.param[1];
        const iz = config.param[2];
        const icat = config.param[2];
        const xCol = dataStore.getRawColumn(ix);
        const yCol = dataStore.getRawColumn(iy);
        const zCol = dataStore.getRawColumn(iz);
        const catCol = dataStore.getRawColumn(icat);
        // seems suboptimal. even if this was better, underlying deck.gl scatterplot data could maybe be Float32Array or something?
        // might consider a custom layer, storying data in texture?
        const indices = new Array(xCol.length).fill().map((_, i) => dataStore.filterArray[i] ? -1 : i).filter(v=> v>-1);
        this.scatterData = new Array(dataStore.filterSize).fill().map((_, j) => {
            const i = indices[j];
            const x = xCol[i];
            const y = yCol[i];
            const z = zCol[i];
            const cat = catCol[i];
            const position = [x, y, z];
            return {
                position, cat, i
            }
        });
        config.scatterData = this.scatterData;
    }
    _update() {
        //const {filterArray} = this.dataStore;
        this.setupScatterplot();
        this.viv.config.scatterData = this.scatterData; //should be redundant (when mutating)
        this.viv._updateProps();
    }
    onDataFiltered(dim) {
        console.log('data filtered', dim);
        super.onDataFiltered(dim);
        this._update();
    }
    onDataAdded(data){
        super.onDataAdded(data);
        this._update();
    }
    onDataChanged(data){
        super.onDataChanged(data);
        this._update();
    }
    onDataHighlighted(data){
        super.onDataHighlighted(data);
        this._update();
    }

    colorByColumn(column) {
        if (!this.scatterData) return;
        this.config.colorby = column;
        const colorFunc = this.getColorFunction(column, true);
        // should allow only the getFillColor function to need change, not the data
        this.viv.config.getScatterFillColor = (d) => colorFunc(d.i);
        // this.setupScatterplot(); // now unnecessary, but still not working without copy:
        /// PJT: WHY OH WHY do I need this in order for it to update props properly?
        this.viv.config.scatterData = [...this.scatterData];
        this.viv._updateProps();
    }
    getColorOptions() {
        if (!this.scatterData) return {};
        return {
            colorby: "all",
            has_default_color: true
        }
    }
    getSettings() {
        return [
            ...super.getSettings(),
            {
                type: "doubleslider",
                label: "clip X",
                min: 0, max: 1,
                current_value: this.viv.clipX,
                continuous: true,
                func: (min, max) => this.viv.setClipX(min, max)
            },
            {
                type: "doubleslider",
                label: "clip Y",
                min: 0, max: 1,
                continuous: true,
                current_value: this.viv.clipY,
                func: (min, max) => this.viv.setClipY(min, max)
            },
            {
                type: "doubleslider",
                label: "clip Z",
                min: 0, max: 1,
                continuous: true,
                current_value: this.viv.clipZ,
                func: (min, max) => this.viv.setClipZ(min, max)
            },
            {
                type: "button",
                label: "Recenter camera",
                func: () => {
                    this.viv.recenterCamera()
                }
            },
        ]
    }
    remove() {
        this.viv.deck?.finalize();
        this.dialog?.close();
        super.remove();
    }
}

const nameOption = {
    type: "string",
    name: "url",
    label: "ome.tiff volume url",
    defaultVal: "/data/t1-head-imj.ome.tiff"
};

const commonOptions = {
    init: (config, dataSource, extraControls) => {
        const {url} = extraControls;
        config.options = {url};
    },
    extra_controls: (dataSource) => {
        return [nameOption];
    }
};

BaseChart.types["viv_volume_scatter_view"] = {
    name: "Viv Volume with Scatterplot",
    class: VivVolume,
    params: [{
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
    {
        type: "text",
        name: "Category Column"
    }
    ],
    ...commonOptions
    // options: [nameOption]
}
BaseChart.types["viv_volume_view"] = {
    name: "Viv Volume View",
    class: VivVolume,
    params: [],
    ...commonOptions
    // options: [nameOption]
}