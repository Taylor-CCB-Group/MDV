import { CompositeLayer, ContourLayer, HeatmapLayer, type Layer, type LayersList } from "deck.gl/typed";
import type { LayerContext, ScatterplotLayerProps } from "deck.gl/typed";
import { ScatterplotExLayer } from './ScatterplotExLayer';
import { Framebuffer } from "@luma.gl/core";
import { ScatterDensityExension } from "./ScatterDeckExtension";
import HeatmapContourExtension, { ExtendableHeatmapLayer } from "./HeatmapContourExtension";


/** In future I think we want something more flexible & expressive,
 * but this should be somewhat compatible with the previous implementation 
 */
export type DualContourLayerProps = {
    contourParameter?: string, //this is param[2] in the original code
    category1?: string,
    category2?: string,
    contour_fill: boolean,
    contour_bandwidth: number,
    contour_intensity: number,
    contour_opacity: number,
}

export type SpatialLayerProps = ScatterplotLayerProps & DualContourLayerProps & {
    // this is how we're doing these category filters for now
    getContourWeight1: (i: number) => number,
    getContourWeight2: (i: number) => number,
};

function rgb(r: number, g: number, b: number, a=255): [number, number, number, number] {
    return [r, g, b, a];
}
// this is not the way to do it...
const contourColors = Array.from({ length: 200 }, (_, i) => {
    const v = i % 20 <= 1 ? 255 : 0;
    return rgb(v, v, v, v);
});
const viridis = [rgb(0, 47, 97), rgb(0, 95, 133), rgb(0, 139, 152), rgb(0, 181, 153), rgb(24, 220, 130), rgb(151, 245, 84), rgb(255, 255, 0)] as const;

export default class SpatialLayer extends CompositeLayer<SpatialLayerProps> {
    static layerName = 'SpatialLayer';
    static defaultProps = ScatterplotExLayer.defaultProps;
    // biome-ignore lint/complexity/noUselessConstructor: tsc error if we remove this?
    constructor(props: SpatialLayerProps) {
        super(props);
    }
    /** "usually not needed for composite layers" but I think this is the right place to setup framebuffers */
    initializeState(context: LayerContext): void {
        super.initializeState(context);
        // this will change with new deck.gl version, context.device rather than context.gl...
        // still considering how to structure coomposite layers...
        // const framebuffer = new Framebuffer(context.gl, {});
    }
    renderLayers(): Layer<SpatialLayerProps> | LayersList {
        // order matters here, we should make a ui where we can easily control it
        return [
            new ScatterplotExLayer(this.getSubLayerProps({
                ...this.props,
                // ...this.props was not including `data`, so we need to add it manually
                // How many more props are going to need manual intervention like this?
                // (could use something like Object.assign())
                // As this develops further, we probably don't want to be using the same props as ScatterplotLayer anyway;
                // sub-layers will be more explicit about what they need.
                data: this.props.data,
                id: 'spatial.scatterplot',
            })),
            // now we need more layers, using gaussian density.
            new ScatterplotExLayer(this.getSubLayerProps({
                ...this.props,
                id: 'spatial.contour1',
                // when deck updates we can use their category filter...
                data: this.props.data,
                opacity: 0.5,
                getRadius: this.props.getContourWeight1,
                updateTriggers: {
                    getRadius: this.props.getContourWeight1
                },
                radiusScale: 100,
                extensions: [new ScatterDensityExension()]
            })),
            // todo applying contour shader, TCM...
            new ExtendableHeatmapLayer(this.getSubLayerProps({
                // ...this.props, // this makes it very slow, we need to be more selective
                getPosition: this.props.getPosition,
                id: 'spatial.heatmap',
                data: this.props.data,
                radiusPixels: 100, //todo - custom version that doesn't have to use pixels (simpler)
                opacity: this.props.contour_intensity,
                getWeight: this.props.getContourWeight2,
                updateTriggers: {
                    getWeight: this.props.getContourWeight2
                },
                // weightsTextureSize: 256,
                colorRange: viridis,
                // colorRange: contourColors,
                // debounceTimeout: 1000,
                // extensions: [new HeatmapContourExtension()] // doesn't do anything here - needs to apply to the right shader
            })),
        ];
    }
}