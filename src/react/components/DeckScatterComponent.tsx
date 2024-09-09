import DeckGL from "@deck.gl/react/typed";
import { observer } from "mobx-react-lite";
import { useChartSize, useFilteredIndices, useParamColumns } from "../hooks";
import { ScatterplotLayer } from "@deck.gl/layers/typed";
import { useMemo, useState } from "react";
import { COORDINATE_SYSTEM } from "deck.gl/typed";
// import { useScatterplotLayer } from "../scatter_state";

export default observer(function DeckScatterComponent() {
    // todo make this non-viv compatible
    // const scatterProps = useScatterplotLayer();
    // const { scatterplotLayer, getTooltip } = scatterProps;
    const [width, height] = useChartSize();
    const [cx, cy] = useParamColumns();
    // const data = useFilteredIndices(); //this fails because of internal assumptions about backgroung_filter etc
    const data = useMemo(() => {
        const data = new Uint32Array(cx.data.length);
        for (let i = 0; i < data.length; i++) {
            data[i] = i;
        }
        return data;
    }, [cx.data]);
    const scatterplotLayer = new ScatterplotLayer({
        id: "scatterplot-layer",
        data,
        pickable: false,
        opacity: 0.8,
        stroked: false,
        filled: true,
        radiusScale: 6,
        radiusMinPixels: 1,
        radiusMaxPixels: 100,
        lineWidthMinPixels: 1,
        getPosition: (i, {target}) => {
            target[0] = cx.data[i];
            target[1] = cy.data[i];
            return target as [number, number];
        },
        getRadius: (d) => 1,
        getFillColor: [0, 140, 240],
        getLineColor: [0, 0, 0],
        coordinateSystem: COORDINATE_SYSTEM.CARTESIAN,
    });
    const [midX, midY, rangeX, rangeY] = useMemo(() => {
        if (data.length === 0) return [0, 0, 1, 1];
        //there is also cx.minMax, cy.minMax - but not for filtered indices
        let minX = Number.POSITIVE_INFINITY;
        let maxX = Number.NEGATIVE_INFINITY;
        let minY = Number.POSITIVE_INFINITY;
        let maxY = Number.NEGATIVE_INFINITY;
        for (let i = 0; i < data.length; i++) {
            const d = data[i];
            const x = cx.data[d];
            const y = cy.data[d];
            if (!Number.isFinite(x) || !Number.isFinite(y)) {
                console.warn("undefined data in scatterplot");
                continue;
            }
            if (x < minX) minX = x;
            if (x > maxX) maxX = x;
            if (y < minY) minY = y;
            if (y > maxY) maxY = y;
        }
        return [(minX + maxX) / 2, (minY + maxY) / 2, maxX - minX, maxY - minY];
    }, [data, cx, cy]);
    const zoom = 4;
    const [viewState, setViewState] = useState<any>({
        width,
        height,
        target: [midX, midY, 0],
        longitude: 0,
        latitude: 0,
        zoom,
        minZoom: -50,
    });

    return (
        <DeckGL 
            layers={[scatterplotLayer]}
            useDevicePixels={true}
            controller={true}
            viewState={viewState}
            onViewStateChange={v => setViewState(v.viewState)}
        />
    );
});