import BaseChart from "../charts/BaseChart.js";
import { createEl } from "../utilities/Elements.js";


class VivVolume extends BaseChart {
    constructor(dataStore, div, config) {
        super(dataStore, div, config, {x:{}, y:{}});
        this.afterAppCreation();
    }
    setSize(x, y) {
        super.setSize(x, y);
        if (this.viv) {
            const b = this._getContentDimensions();
            this.viv.setSize(b.width, b.height, this.app);
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
        import('../webgl/VivViewer.js').then(({default:VivViewer}) => {
            const vivConfig = {
                url: '/data/t1-head-imj.ome.tiff',
                use3d: true
            }
            // not really what we want...
            const iv = {
                x_scale: 1, y_scale: 1, offset: 0
            }
            this.viv = new VivViewer(this.vivCanvas, vivConfig, iv);
        });
        // return super.afterAppCreation();
    }
}

BaseChart.types["viv_volume_view"] = {
    name: "Viv Volume View",
    class: VivVolume,
    params:[
        // hmmm, url isn't a column...
        // {
        //     type: "string",
        //     name: "url"
        // }
    ]
}