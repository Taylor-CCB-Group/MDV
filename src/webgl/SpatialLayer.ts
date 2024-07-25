import { CompositeLayer, type Layer, type LayersList } from "deck.gl/typed";
import type { LayerContext, ScatterplotLayerProps } from "deck.gl/typed";
import { ScatterplotExLayer } from './ScatterplotExLayer';
import { Framebuffer } from "@luma.gl/core";
import { ScatterDensityExension } from "./ScatterDeckExtension";


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
        const framebuffer = new Framebuffer(context.gl, {});
    }
    renderLayers(): Layer<SpatialLayerProps> | LayersList {
        // for this to be rendered in a VivViewer, the id must have an appropriate pattern
        // to satisfy `layer.id.includes(getVivId(viewport.id));` in VivViewer.tsx (we can use MDVivViewer for testing)
        const id = `spatial_composite.${this.props.id}`;
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
            // todo rendering into intermediate buffers, applying contour shader, blending them together...
        ];
    }
}