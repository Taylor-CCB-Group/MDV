import DeckGL from "@deck.gl/react/typed";
import { OrthographicView } from '@deck.gl/core';
import { observer } from "mobx-react-lite";
import { useChartSize, useConfig, useFilteredIndices, useParamColumns } from "../hooks";
import { ScatterplotLayer } from "@deck.gl/layers/typed";
import { useEffect, useId, useMemo, useState } from "react";
import type { ScatterPlotConfig } from "../scatter_state";
import { useChart } from "../context";
// import { useScatterplotLayer } from "../scatter_state";

export default observer(function DeckScatterComponent() {
    const id = useId();
    const [width, height] = useChartSize();
    const [cx, cy] = useParamColumns();
    const data = useFilteredIndices();
    const config = useConfig<ScatterPlotConfig>();
    const { opacity, radius, course_radius } = config;
    const radiusScale = radius * course_radius;
    const chart = useChart();
    const colorBy = (chart as any).colorBy;


    const scatterplotLayer = new ScatterplotLayer({
        id: `scatterplot-layer-${id}`,
        data,
        pickable: true,
        opacity,
        stroked: false,
        filled: true,
        radiusScale,
        getPosition: (index, {target}) => {
            target[0] = cx.data[index];
            target[1] = cy.data[index];
            return target as [number, number];
        },
        getFillColor: colorBy ?? [255, 255, 255],
        getLineColor: [0, 0, 0],
        updateTriggers: {
            getFillColor: colorBy,
        }
    });
    const [viewState, setViewState] = useState<any>({
        width,
        height,
        target: [0, 0, 0],
        zoom: 0,
        minZoom: -50,
    });

    //this should be factored out into a hook
    useEffect(() => {
        if (data.length === 0) return;// [0, 0, 1, 1];
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
            width,
            height,
            target: [(minX + maxX) / 2, (minY + maxY) / 2, 0],
            zoom: Math.log2(Math.min(width/(maxX - minX), height/(maxY - minY))) - 0.1,
        });
        //return [(minX + maxX) / 2, (minY + maxY) / 2, maxX - minX, maxY - minY];
    }, [data, cx, cy, width, height]);

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