import { CompositeLayer, type Layer, type LayerProps, type LayersList } from "@deck.gl/core";
import type { ScatterplotLayerProps } from "@deck.gl/layers";
import { ScatterplotLayer } from "@deck.gl/layers";
import { DensityLayer } from "@intraspatial/deck";
import { HeatmapLayer } from "deck.gl";
import type { ContourLayerProps } from "@/react/contour_state";
import { TriangleLayerContours } from "./HeatmapContourExtension";
import { contourToDensityProps } from "./intraspatialDensity";
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

/**
 * The density sublayers for `contourLayers`. `subLayerProps` is the parent composite's
 * `getSubLayerProps`, which supplies the model matrix and the other props deck forwards.
 *
 * TriangleLayerContours is GLSL-only. On WebGPU, contours are @intraspatial/deck's
 * DensityLayer instead: KDE splatting in data units with a WGSL presentation shader,
 * converted from the heatmap props in `intraspatialDensity.ts`.
 */
function densitySubLayers(
    contourLayers: ContourLayerProps[],
    isWebGPU: boolean,
    // biome-ignore lint/suspicious/noExplicitAny: deck types getSubLayerProps' result as any.
    subLayerProps: (props: ContourLayerProps) => any,
): Layer[] {
    const contourSubLayerProps = {
        triangle: {
            type: TriangleLayerContours,
        },
    };
    return contourLayers
        .filter((l) => l)
        .map((props) => {
            const { extensions, ...p } = subLayerProps(props);
            if (isWebGPU) return new DensityLayer(contourToDensityProps(p));
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
}

export type DensityOverlayLayerProps = Pick<LayerProps, "id" | "modelMatrix" | "pickable"> & {
    contourLayers: ContourLayerProps[];
};

/**
 * Density contours on their own, for a chart that lists them as a layer apart from the
 * points (the SpatialData chart's "Density"). Other charts draw them inside `SpatialLayer`.
 */
export class DensityOverlayLayer extends CompositeLayer<DensityOverlayLayerProps> {
    static layerName = "DensityOverlayLayer";
    renderLayers(): LayersList {
        return densitySubLayers(this.props.contourLayers, isWebGPULayer(this), (p) => this.getSubLayerProps(p));
    }
}

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
        const densityLayers = densitySubLayers(contourLayers, isWebGPULayer(this), (p) => this.getSubLayerProps(p));
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
