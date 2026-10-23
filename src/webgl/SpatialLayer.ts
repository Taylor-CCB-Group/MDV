import { CompositeLayer, type Layer, type LayersList } from "@deck.gl/core";
import type { ScatterplotLayerProps } from "@deck.gl/layers";
import { ScatterplotLayer } from "@deck.gl/layers";
import { TriangleLayerContours } from "./HeatmapContourExtension";
import type { ContourLayerProps } from "@/react/contour_state";
import { HeatmapLayer } from "deck.gl";
import { MdvScatterplotLayer } from "./MdvScatterplotLayer";
import { isWebGPULayer } from "./wgslLayerExtension";

export type SpatialLayerProps = ScatterplotLayerProps & {
    //pending typing that allows for other kinds of layers etc
    contourLayers: ContourLayerProps[];
    highlightedData?: number[];
    highlightLineWidth?: number;
    highlightRadiusScale?: number;
    highlightLineColor?: [number, number, number];
};

export default class SpatialLayer extends CompositeLayer<SpatialLayerProps> {
    static layerName = "SpatialLayer";
    static defaultProps = ScatterplotLayer.defaultProps;
    // biome-ignore lint/complexity/noUselessConstructor: tsc error if we remove this?
    constructor(props: SpatialLayerProps) {
        super(props);
    }
    renderLayers(): Layer<SpatialLayerProps> | LayersList {
        // order matters here, we should make a ui where we can easily control it
        const {
            contourLayers,
            highlightedData = [],
            highlightLineWidth = 1,
            highlightRadiusScale,
            highlightLineColor = [255, 255, 255],
        } = this.props;
        // TriangleLayerContours is GLSL-only. On WebGPU the heatmap keeps deck's own
        // (WGSL) triangle layer: plain heatmap shading, no contour lines. Contours are
        // expected to move to KDE splatting with a separate presentation shader rather
        // than be ported as they are.
        const contourSubLayerProps = isWebGPULayer(this)
            ? {}
            : {
                  triangle: {
                      type: TriangleLayerContours,
                  },
              };
        const densityLayers = contourLayers.filter((l) => l).map(props => {
            const { extensions, ...p } = this.getSubLayerProps(props);
            return new HeatmapLayer({
                ...p,
                _subLayerProps: {
                    ...contourSubLayerProps,
                    "triangle-layer": {
                        contourOpacity: p.contourOpacity,
                        contourFill: p.contourFill,
                        fillOpacity: p.fillOpacity,
                    },
                },
            });
        });
        return [
            // add 'grey-out' layer here... that implies a different type of data being passed.
            // it might want to know about background_filter...

            new MdvScatterplotLayer(
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
            highlightedData.length > 0
                ? new MdvScatterplotLayer(
                      this.getSubLayerProps({
                          ...this.props,
                          data: highlightedData,
                          id: "spatial.highlight",
                          filled: false,
                          stroked: true,
                          pickable: false,
                          radiusScale: highlightRadiusScale ?? this.props.radiusScale,
                          getLineWidth: highlightLineWidth,
                          getLineColor: highlightLineColor,
                          lineWidthMinPixels: 1,
                      }),
                  )
                : null,
            // consider trying to render density map at lower resolution, then upscaling on debounce.
            ...densityLayers,
        ].filter(Boolean) as LayersList;
    }
}
