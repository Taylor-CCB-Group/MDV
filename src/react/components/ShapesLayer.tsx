import { GeoJsonLayer } from "@deck.gl/layers";
import { DataFilterExtension } from "@deck.gl/extensions";
import { useMemo } from "react";
import type { Feature, GeoJsonProperties, Geometry } from "geojson";
import { useProject } from "@/modules/ProjectContext";
import { useChartID, useFieldSpec, useFilteredIndices, useRegion } from "../hooks";
import { useHighlightedIndices, useHighlightRows } from "../selectionHooks";
import { getVivId } from "./avivatorish/MDVivViewer";

type GeoJsonFeatureProperties = GeoJsonProperties & {
    id?: string | number | null;
    DN?: string | number | null;
};

type GeoJsonFeature = Feature<Geometry, GeoJsonFeatureProperties>;

function getGeoJsonFeatureId(feature: GeoJsonFeature | undefined): string | null {
    const featureId = feature?.id ?? feature?.properties?.id;
    if (featureId == null) return null;
    return String(featureId);
}

export function useShapesLayer(showShapes: boolean) {
    const id = useChartID();
    const { root } = useProject();
    const { json } = useRegion();
    const cellIdColumn = useFieldSpec("EntityID"); //!!!!!! we need a better basis on which to associate id. cell_id has `_N` after merge...
    const filteredIndices = useFilteredIndices();
    const highlightedIndices = useHighlightedIndices();
    const highlightRows = useHighlightRows();
    const layerId = `json_${getVivId(`${id}detail-react`)}`;
    const dataUrl = useMemo(() => (json ? `${root}/${json}` : null), [json, root]);

    const geoJsonState = useMemo(() => {
        const cellIds = cellIdColumn?.data;
        if (!cellIds) return null;
        const visibleIds = new Set<string>();
        const highlightedIds = new Set<string>();
        const rowIndexById = new Map<string, number>();

        for (let rowIndex = 0; rowIndex < cellIds.length; rowIndex += 1) {
            const cellId = cellIdColumn.getValue(rowIndex);
            if (cellId == null) continue;
            const key = String(cellId);
            if (!rowIndexById.has(key)) rowIndexById.set(key, rowIndex);
        }
        for (const rowIndex of filteredIndices) {
            const cellId = cellIdColumn.getValue(rowIndex);
            if (cellId != null) {
                visibleIds.add(String(cellId));
            }
        }
        for (const rowIndex of highlightedIndices) {
            const cellId = cellIdColumn.getValue(rowIndex);
            if (cellId != null) highlightedIds.add(String(cellId));
        }

        return { visibleIds, highlightedIds, rowIndexById };
    }, [cellIdColumn?.data, filteredIndices, highlightedIndices]);

    return useMemo(() => {
        type GeoJsonLayerProps = ConstructorParameters<typeof GeoJsonLayer<GeoJsonFeature>>[0];
        if (!dataUrl) return null;

        const layerProps: GeoJsonLayerProps & {
            getFilterValue: (feature: GeoJsonFeature) => number;
            filterRange: [number, number];
        } = {
            id: layerId,
            data: dataUrl,
            opacity: 1,
            filled: true,
            getFillColor: (feature: GeoJsonFeature) => {
                const featureId = getGeoJsonFeatureId(feature);
                if (featureId != null && geoJsonState?.highlightedIds.has(featureId)) {
                    return [255, 191, 0, 180];
                }
                return [255, 255, 255, 110];
            },
            getLineColor: (feature: GeoJsonFeature) => {
                const featureId = getGeoJsonFeatureId(feature);
                if (featureId != null && geoJsonState?.highlightedIds.has(featureId)) {
                    return [255, 191, 0, 255];
                }
                return [255, 255, 255, 170];
            },
            getLineWidth: (feature: GeoJsonFeature) => {
                const featureId = getGeoJsonFeatureId(feature);
                if (featureId != null && geoJsonState?.highlightedIds.has(featureId)) return 3;
                return 2;
            },
            lineWidthMinPixels: 1,
            getPointRadius: 10,
            pickable: true,
            autoHighlight: true,
            extensions: [new DataFilterExtension({ filterSize: 1 })],
            getFilterValue: (feature: GeoJsonFeature) => {
                const featureId = getGeoJsonFeatureId(feature);
                if (featureId == null || !geoJsonState) return 1;
                return geoJsonState.visibleIds.has(featureId) ? 1 : 0;
            },
            filterRange: [0.5, 1],
            updateTriggers: {
                getFillColor: [geoJsonState?.highlightedIds],
                getLineColor: [geoJsonState?.highlightedIds],
                getLineWidth: [geoJsonState?.highlightedIds],
                getFilterValue: [geoJsonState?.visibleIds, geoJsonState?.highlightedIds],
            },
            onClick: ({ object }) => {
                const featureId = getGeoJsonFeatureId(object as GeoJsonFeature | undefined);
                if (featureId == null) return;
                const rowIndex = geoJsonState?.rowIndexById.get(featureId);
                if (rowIndex === undefined) return;
                highlightRows([rowIndex]);
            },
            //@ts-expect-error GeoJson getText: might think about using zod to type/validate this
            getText: (feature) => feature.properties?.DN,
            getTextColor: [255, 255, 255, 255],
            getTextSize: 12,
            textBackground: true,
            visible: showShapes,
        };
        return new GeoJsonLayer<GeoJsonFeature>(layerProps);
    }, [dataUrl, geoJsonState, highlightRows, layerId, showShapes]);
}

export default useShapesLayer;
