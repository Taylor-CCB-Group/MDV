import { WGL2DI } from "../webgl/WGL2DI.js";
import WGLChart from "./WGLChart.js";
import BaseChart from "./BaseChart.js";
import { createEl } from "../utilities/Elements.js";
import { ImageArray } from "../webgl/ImageArray.js";


class ImageScatterChart extends BaseChart {
    app: WGL2DI;
    canvas: HTMLCanvasElement;
    constructor(dataStore, div, config) {
        super(dataStore, div, config);
        const canvas = this.canvas = createEl("canvas", {}, this.contentDiv);
        //const gl = canvas.getContext("webgl2"); // do we need to take care of disposing resources as well?
        const imageArray = new ImageArray(dataStore, canvas, {
            base_url: config.images.base_url,
            image_type: "png",
            image_key: config.images.image_key, //"figure_id",
            width: 512, //TODO: specify this properly
            height: 512,
        });

        // when we implement for real, we'll need to do something like this
        // with a new 'mode' parameter...
        // this.app = new WGL2DI(div, config);
    }
}

BaseChart.types["ImageScatterChart"] = {
    class: ImageScatterChart,
    name: "Image Scatter Chart",
    required: ["images"],
    

    init: (config, dataSource, extraControls) => {
        //get the available images
        const i = dataSource.images[extraControls.image_set];
        console.log('ImageScatterChart param', config.param);
        //set the base url and type
        config.images = {
            base_url: i.base_url,
            type: "png", //todo: allow this to be specified
            image_key: i.key_column, //nb ImageTableChart has this as config.param
        }
        config.sortBy = extraControls.sort_by;
        config.sortOrder = extraControls.sort_order;
    },
    extra_controls: (dataSource) => {
        const imageSets = [];
        for (let iname in dataSource.images) {
            imageSets.push({ name: iname, value: iname })
        }
        console.log('imageSets', imageSets);
        const sortableColumns = dataSource.getLoadedColumns().map(c => ({ name: c, value: c }));
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
    ]
};

export default ImageScatterChart;