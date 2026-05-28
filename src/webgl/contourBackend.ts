import type { Layer } from "@deck.gl/core";
import type { Device } from "@luma.gl/core";
import { HeatmapLayer } from "deck.gl";
import { TriangleLayerContours } from "./HeatmapContourExtension";
import type { ContourLayerProps } from "@/react/contour_state";

export type ContourBackendStaticResources = {
    device: Device;
    width: number;
    height: number;
};

export type ContourBackend = {
    id: string;
    buildStaticResources?: (
        resources: ContourBackendStaticResources,
    ) => Record<string, unknown> | null;
    buildDynamicLayer: (args: {
        contourLayer: NonNullable<ContourLayerProps>;
        id: string;
        radiusPixels: number;
        debounce?: number;
        weightsTextureSize?: number;
    }) => Layer;
    dispose?: (resources?: Record<string, unknown> | null) => void;
};

function buildHeatmapContourLayer(
    props: {
        id: string;
        radiusPixels: number;
        debounce?: number;
        weightsTextureSize?: number;
        viewId?: string;
    } & Record<string, unknown>,
): Layer {
    const { extensions, ...p } = props as Record<string, unknown>;
    return new HeatmapLayer({
        ...(p as object),
        _subLayerProps: {
            triangle: {
                type: TriangleLayerContours,
            },
            "triangle-layer": {
                contourOpacity: p.contourOpacity as number,
                contourFill: p.contourFill as number,
                fillOpacity: p.fillOpacity as number,
            },
        },
    }) as unknown as Layer;
}

export const heatmapContourBackend: ContourBackend = {
    id: "heatmap-contour",
    buildDynamicLayer: ({
        contourLayer,
        id,
        radiusPixels,
        debounce = 250,
        weightsTextureSize = 256,
    }) =>
        buildHeatmapContourLayer({
            ...contourLayer,
            id,
            radiusPixels,
            debounce,
            weightsTextureSize,
        }),
};

export const framebufferContourBackendPlaceholder: ContourBackend = {
    id: "framebuffer-contour-placeholder",
    buildDynamicLayer: ({ contourLayer, id, radiusPixels }) =>
        buildHeatmapContourLayer({
            ...contourLayer,
            id,
            radiusPixels,
            debounce: 250,
            weightsTextureSize: 256,
        }),
};
