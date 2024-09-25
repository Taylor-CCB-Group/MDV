import { CompositeLayer, type Layer, type LayersList } from "deck.gl/typed";
import type { LayerContext, ScatterplotLayerProps } from "deck.gl/typed";
import { ScatterplotExLayer } from "./ScatterplotExLayer";
import { Framebuffer } from "@luma.gl/core";
import { ScatterDensityExension } from "./ScatterDeckExtension";
import HeatmapContourExtension, {
    ExtendableHeatmapLayer,
} from "./HeatmapContourExtension";
import type { useContour } from "@/react/contour_state";

export type SpatialLayerProps = ScatterplotLayerProps & {
    //pending typing that allows for other kinds of layers etc
    contourLayers: ReturnType<typeof useContour>[];
};

export default class SpatialLayer extends CompositeLayer<SpatialLayerProps> {
    static layerName = "SpatialLayer";
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
            new ScatterplotExLayer(
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
            // now we need more layers, using gaussian density.
            // consider trying to render density map at lower resolution, then upscaling on debounce.
            ...this.props.contourLayers
                .filter((l) => l)
                .map(
                    (props) =>
                        new ExtendableHeatmapLayer(
                            this.getSubLayerProps(props),
                        ),
                ),
        ];
    }
}
