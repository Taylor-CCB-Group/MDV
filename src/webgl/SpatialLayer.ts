import { CompositeLayer, type Layer, type LayersList } from "deck.gl/typed";
import type { ScatterplotLayerProps } from "deck.gl/typed";
import { ScatterplotExLayer } from './ScatterplotExLayer';
import { getVivId } from "@/react/components/avivatorish/MDVivViewer";


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
    renderLayers(): Layer<SpatialLayerProps> | LayersList {
        // for this to be rendered in a VivViewer, the id must have an appropriate pattern
        // to satisfy `layer.id.includes(getVivId(viewport.id));` in VivViewer.tsx (we can use MDVivViewer for testing)
        const id = `spatial_composite.${this.props.id}`;
        return [
            new ScatterplotExLayer({
                ...this.props,
                // I don't understand why this is necessary, but was ending up with data being [] down the line.
                // How many more props are going to need manual intervention like this?
                // As this develops further, we probably don't want to be using the same props as ScatterplotLayer anyway;
                // sub-layers will be more explicit about what they need.
                data: this.props.data,
                id,
            }),
        ];
    }
}