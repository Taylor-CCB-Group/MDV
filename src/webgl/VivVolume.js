import BaseChart from "../charts/BaseChart.js";
import { ColorChannelDialog } from "../charts/VivScatterPlot.js";
import { createEl } from "../utilities/Elements.js";
import VivViewer from '../webgl/VivViewer.js';

class VivVolume extends BaseChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config, {x:{}, y:{}});
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
            url: this.config.options.url, //'/data/t1-head-imj.ome.tiff',
            use3d: true
        }
        // not really what we want... basically ignored with use3d.
        const iv = {
            x_scale: 1, y_scale: 1, offset: 0
        }
        // we may not necessarily want to re-use VivViewer?
        // although some common features with channels etc.
        this.viv = new VivViewer(this.vivCanvas, vivConfig, iv);
    }

    getSettings(conf) {
        return [
            {
                type: "button",
                label: "Recenter camera",
                func: () => {
                    this.viv.recenterCamera()
                }
            }
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
    options: [nameOption]
}
BaseChart.types["viv_volume_view"] = {
    name: "Viv Volume View",
    class: VivVolume,
    params: [],
    options: [nameOption]
}