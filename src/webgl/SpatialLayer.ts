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
    contourParameter: string,
    category1?: string,
    category2?: string,
    contour_bandwidth: number,
    contour_intensity: number,
    contour_opacity: number,
}

export type SpatialLayerProps = ScatterplotLayerProps & DualContourLayerProps;

function rgb(r: number, g: number, b: number, a=255): [number, number, number, number] {
    return [r, g, b, a];
}
// this is not the way to do it...
const contourColors = Array.from({ length: 200 }, (_, i) => {
    const v = i % 20 <= 1 ? 255 : 0;
    return rgb(v, v, v, v);
});
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
        // const framebuffer = new Framebuffer(context.gl, {});
    }
    renderLayers(): Layer<SpatialLayerProps> | LayersList {
        // for this to be rendered in a VivViewer, the id must have an appropriate pattern
        // to satisfy `layer.id.includes(getVivId(viewport.id));` in VivViewer.tsx (we can use MDVivViewer for testing)
        const id = `spatial_composite.${this.props.id}`;
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
                opacity: this.props.contour_opacity,
                radiusScale: 100,
                extensions: [new ScatterDensityExension()]
            })),
            // todo applying contour shader, TCM...
            new ExtendableHeatmapLayer(this.getSubLayerProps({
                ...this.props,
                id: 'spatial.heatmap',
                data: this.props.data,
                radiusPixels: 100, //todo - custom version that doesn't have to use pixels (simpler)
                opacity: 0.2,
                // weightsTextureSize: 256,
                colorRange: [rgb(0, 47, 97), rgb(0, 95, 133), rgb(0, 139, 152), rgb(0, 181, 153), rgb(24, 220, 130), rgb(151, 245, 84), rgb(255, 255, 0)],
                // colorRange: contourColors,
                debounceTimeout: 1000,
                extensions: [new HeatmapContourExtension()]
            })),
            // nb, ContourLayer exists, but I think I want to do something different...
            // that said, I should at least see how it works...
            // new ContourLayer(this.getSubLayerProps({
            //     ...this.props,
            //     id: 'spatial.contour2',
            //     data: this.props.data,
            //     // opacity: this.props.contour_opacity,
            //     // radiusScale: 100, // I don't understand how to embiggen the influence of the points
            //     cellSize: 1e6, //<<< very important, it can't infer sensible defaults
            //     // getWeight: (d) => 10,
            //     contours: [
            //         { threshold: 0.1, color: [255, 0, 0] },
            //         // { threshold: 0.5, color: [0, 255, 0] },
            //         { threshold: [0.5, 1.5], color: [255, 255, 255] },
            //         { threshold: 2, color: [255, 255, 255] },
            //     ],
            //     opacity: 0.1,
            // })),
        ];
    }
}