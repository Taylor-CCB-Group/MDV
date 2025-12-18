import { CompositeLayer, Layer, type LayersList } from "@deck.gl/core";
import type { LayerContext } from "@deck.gl/core";
import { type ScatterplotLayerProps, ScatterplotLayer } from "@deck.gl/layers";
import { TriangleLayerContours } from "./HeatmapContourExtension";
import type { ContourLayerProps } from "@/react/contour_state";
import { HeatmapLayer } from "deck.gl";
import type { Framebuffer } from "@luma.gl/core";

export type SpatialLayerProps = ScatterplotLayerProps & {
    //pending typing that allows for other kinds of layers etc
    contourLayers: ContourLayerProps[];
};
type DensityLayerProps = ScatterplotLayerProps & {
    // are framebuffers appropriate for use as props?
    framebuffer: Framebuffer;
};
class DensityLayer extends Layer<DensityLayerProps> {
    layerName = "DensityLayer";
    static defaultProps = {};
    constructor(props: DensityLayerProps) {
        super(props);
        console.log("DensityLayer constructor, very much under construction");
    }
    initializeState(context: LayerContext) {
        const { device } = context;
        // how should we tell what dimensions to use?
        const size = {
            width: 256,
            height: 256,
        };
        //https://luma.gl/docs/api-reference/core/resources/framebuffer/
        //maybe we can have a composite layer with a '2d-array' that gets passed in...
        console.log(`Initialising DensityLayer '${this.props.id}'`);
        const framebuffer = device.createFramebuffer({
            ...size,
            colorAttachments: [
                device.createTexture({ format: "r32float", ...size }),
            ],
        });
        this.setState({ framebuffer });
    }
    draw() {
        // const { framebuffer } = this.props;
        const models = this.getModels();
        // models.forEach(model => model.);
    }
}

export default class SpatialLayer extends CompositeLayer<SpatialLayerProps> {
    static layerName = "SpatialLayer";
    static defaultProps = ScatterplotLayer.defaultProps;
    static bufferExperiment = false;
    // biome-ignore lint/complexity/noUselessConstructor: tsc error if we remove this?
    constructor(props: SpatialLayerProps) {
        super(props);
    }
    /** "usually not needed for composite layers" but I think this is the right place to setup framebuffers */
    initializeState(context: LayerContext): void {
        super.initializeState(context);
        if (!SpatialLayer.bufferExperiment) return;
        // this will change with new deck.gl version, context.device rather than context.gl...
        // still considering how to structure coomposite layers...
        // const framebuffer = new Framebuffer(context.gl, {});
        // many ways to skin a cat here... want to render each contourLayer into a lower resolution (float?) fbo,
        // with a gaussian scatterplot to be used as a density map
        // that should then be fed into a shader that can take these densities and do some calculations
        // const framebuffer = this.context.device.createFramebuffer({width: 256, height: 256})
        const { device } = context;
        const size = {
            width: 256, //tbd, some fraction of the original size
            height: 256,
            depth: 256, //<- should be the number of fields/layers
        };
        const texture = device.createTexture({
            // dimension: "3d",
            format: "r32float",
            ...size,
        });
        const framebuffer = device.createFramebuffer({
            colorAttachments: [texture],
        });
    }
    renderLayers(): Layer<SpatialLayerProps> | LayersList {
        // order matters here, we should make a ui where we can easily control it
        const { contourLayers } = this.props;
        // future work:
        // const densityLayers = contourLayers.map((layer) => {
        //     const { extensions, ...p } = this.getSubLayerProps(layer);
        //     return new DensityLayer({
        //         ...p,
        //     });
        // });
        return [
            // add 'grey-out' layer here... that implies a different type of data being passed.
            // it might want to know about background_filter...

            // now we need more layers, using gaussian density.
            // consider trying to render density map at lower resolution, then upscaling on debounce.
            ...contourLayers
                .filter((l) => l)
                .map((props) => {
                    const { extensions, ...p } = this.getSubLayerProps(props);
                    // todo: maybe encapsulate this in a different way
                    //patch so that it doesn't try to use incompatible extensions used by the ScatterplotLayer
                    //up for review...
                    if (extensions && extensions.length > 0) {
                        console.log(
                            "pending review how extensions interact with sublayers - filtering out from subLayerProps",
                            extensions,
                        );
                    }
                    return new HeatmapLayer({
                        ...p,
                        // extensions: [],
                        _subLayerProps: {
                            triangle: {
                                type: TriangleLayerContours,
                            },
                            "triangle-layer": {
                                contourOpacity: p.contourOpacity,
                            },
                        },
                    });
                }),
            new ScatterplotLayer(
                this.getSubLayerProps({
                    ...this.props,
                    // ...this.props was not including `data`, so we need to add it manually
                    // How many more props are going to need manual intervention like this?
                    // (could use something like Object.assign())
                    // As this develops further, we probably don't want to be using the same props as ScatterplotLayer anyway;
                    // sub-layers will be more explicit about what they need.
                    data: this.props.data,
                    id: "spatial.scatterplot",
                }),
            ),
            // ...densityLayers,
        ];
    }
    draw(opts: any) {
        super.draw(opts);
        //? draw sublayers into framebuffer(s) and then use that in rendering with a custom shader
    }
}
