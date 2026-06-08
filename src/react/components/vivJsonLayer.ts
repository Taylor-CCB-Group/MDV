import { GeoJsonLayer } from "@deck.gl/layers";
import { useMemo } from "react";
import { useChartID, useRegion } from "../hooks";
import { useProject } from "@/modules/ProjectContext";
import { getVivId } from "./avivatorish/MDVivViewer";
import { tagDeckLayerViewportScope } from "./deckLayerViewportScope";

/** Overlay uses `detail-react` view id; chart-array grid uses a distinct deck layer id after overlay unmount. */
export type JsonLayerDeckContext = "overlay" | "chart-array";

export function useJsonLayer(
    showJson: boolean,
    deckContext: JsonLayerDeckContext = "overlay",
) {
    const id = useChartID();
    const { root } = useProject();
    const { json } = useRegion(); // return type is 'any' and we assume 'json' will be a string - but want that to be different in future.
    const idSuffix = deckContext === "chart-array" ? "chart-array-grid" : "detail-react";
    const layer_id = `json_${getVivId(`${id}${idSuffix}`)}`;
    const layer = useMemo(() => {
        if (!json) return null;
        const geoJsonLayer = new GeoJsonLayer({
            id: layer_id,
            data: `${root}/${json}`,
            opacity: 1,
            filled: true,
            getFillColor: (f) => [255, 255, 255, 150],
            getLineColor: (f) => [255, 255, 255, 150],
            getLineWidth: 2,
            lineWidthMinPixels: 1,
            getPointRadius: 10,
            pickable: true,
            autoHighlight: true,
            //@ts-expect-error GeoJson getText: might think about using zod to type/validate this
            getText: (f) => f.properties.DN,
            getTextColor: [255, 255, 255, 255],
            getTextSize: 12,
            textBackground: true,
            visible: showJson,
        });
        return tagDeckLayerViewportScope(geoJsonLayer, "chart-shared");
    }, [json, showJson, layer_id, root]);
    return layer;
}
