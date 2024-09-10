import DeckGL from "@deck.gl/react/typed";
import { OrthographicView } from '@deck.gl/core';
import { observer } from "mobx-react-lite";
import { useChartSize, useFilteredIndices, useParamColumns } from "../hooks";
import { ScatterplotLayer } from "@deck.gl/layers/typed";
import { useMemo, useState } from "react";
import { COORDINATE_SYSTEM } from "deck.gl/typed";
// import { useScatterplotLayer } from "../scatter_state";

export default observer(function DeckScatterComponent() {
    const [width, height] = useChartSize();
    const [cx, cy] = useParamColumns();
    const data = useFilteredIndices();
    // const data = useMemo(() => ({ length: cx.data.length }), [cx.data.length]);
    const scatterplotLayer = new ScatterplotLayer({
        id: "scatterplot-layer",
        data,
        pickable: false,
        opacity: 0.8,
        stroked: false,
        filled: true,
        radiusScale: 10,
        radiusMinPixels: 0.1,
        radiusMaxPixels: 100,
        lineWidthMinPixels: 1,
        getPosition: (index, {target}) => {
            target[0] = cx.data[index];
            target[1] = cy.data[index];
            return target as [number, number];
        },
        getRadius: (d) => 1,
        getFillColor: [0, 140, 240],
        getLineColor: [0, 0, 0],
        // coordinateSystem: COORDINATE_SYSTEM.CARTESIAN,
    });
    const [viewState, setViewState] = useState<any>({
        width,
        height,
        target: [0, 0, 0],
        zoom: 0,
        minZoom: -50,
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
        setViewState({
            ...viewState,
            target: [(minX + maxX) / 2, (minY + maxY) / 2, 0],
            zoom: Math.log2(Math.min(width/(maxX - minX), height/(maxY - minY))) - 0.1,
        });
        return [(minX + maxX) / 2, (minY + maxY) / 2, maxX - minX, maxY - minY];
    }, [data, cx, cy]);

    // we need an OrthographicView to prevent wrapping etc...
    const view = useMemo(() => new OrthographicView({
        id: "scatterplot-view",
        controller: true,
        width,
        height,
        x: 0,
        y: 0,
    }), [width, height]);

    return (
        <>
        <DeckGL 
            layers={[scatterplotLayer]}
            useDevicePixels={true}
            controller={true}
            viewState={viewState}
            views={view}
            onViewStateChange={v => setViewState(v.viewState)}
        />
        </>
    );
});