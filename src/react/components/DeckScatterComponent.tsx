import DeckGL from "@deck.gl/react";
import { OrthographicView, OrbitView } from '@deck.gl/core';
import { observer } from "mobx-react-lite";
import { useChartSize, useConfig, useFilteredIndices, useParamColumns } from "../hooks";
import { ScatterplotLayer } from "@deck.gl/layers";
import { useEffect, useId, useMemo } from "react";
import { useChart } from "../context";
import type { DeckScatterConfig } from "./DeckScatterReactWrapper";
import { action } from "mobx";
import type { OrbitViewState } from "@deck.gl/core";


/** todo this should be common for viv / scatter_state, pending refactor */
function useZoomOnFilter() {
    const [width, height] = useChartSize();
    const [cx, cy, cz] = useParamColumns();
    const data = useFilteredIndices();
    const config = useConfig<DeckScatterConfig>();
    const chart = useChart() as any;
    const { pendingRecenter } = chart;
    useEffect(() => {
        if (chart.ignoreStateUpdate) return;
        if (data.length === 0) return;// [0, 0, 1, 1];
        if (!pendingRecenter && !config.zoom_on_filter) return;
        //there is also cx.minMax, cy.minMax - but not for filtered indices
        let minX = Number.POSITIVE_INFINITY;
        let maxX = Number.NEGATIVE_INFINITY;
        let minY = Number.POSITIVE_INFINITY;
        let maxY = Number.NEGATIVE_INFINITY;
        // todo z... probably refactor for less repetition
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
        // setViewState(
        action(() => {
            chart.ignoreStateUpdate = true;
            chart.pendingRecenter = false;
            const oldViewState = config.viewState;
            config.viewState = {
                target: [(minX + maxX) / 2, (minY + maxY) / 2, 0],
                zoom: Math.log2(Math.min(width / (maxX - minX), height / (maxY - minY))) - 0.1,
            };
            if (config.dimension === "3d") {
                // a bit clunktastic, might make more fancy typescript-y way to do this
                // (thanks copilot for the word clunktastic - maybe a little strong)
                const vs = oldViewState as OrbitViewState;
                const newVs = config.viewState as OrbitViewState;
                if (vs) {
                    newVs.rotationOrbit = vs.rotationOrbit || 0;
                    newVs.rotationX = vs.rotationX || 0;
                }
            }

            chart.ignoreStateUpdate = false;
        })();
        //return [(minX + maxX) / 2, (minY + maxY) / 2, maxX - minX, maxY - minY];
    }, [data, cx, cy, width, height, config, pendingRecenter, chart]);
}



export default observer(function DeckScatterComponent() {
    const id = useId();
    const [width, height] = useChartSize();
    const [cx, cy, cz] = useParamColumns();
    const data = useFilteredIndices();
    const config = useConfig<DeckScatterConfig>();
    const { opacity, radius, course_radius, viewState } = config;
    // const is2d = dimension === "2d";
    //todo more clarity on radius units - but large radius was causing big problems after deck upgrade
    const radiusScale = radius * course_radius * 0.01;// * (is2d ? 1 : 0.01);
    const chart = useChart();
    const colorBy = (chart as any).colorBy;

    useZoomOnFilter();


    const scatterplotLayer = useMemo(() => new ScatterplotLayer({
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
            if (cz) target[2] = cz.data[index];
            return target as [number, number];
        },
        getFillColor: colorBy ?? [61, 126, 180],
        getLineColor: [0, 0, 0],
        updateTriggers: {
            getFillColor: colorBy,
        },
        billboard: true,
        parameters: {
            depthTest: false,
        }
    }), [data, cx, cy, cz, colorBy, opacity, radiusScale, id]);

    // we need an OrthographicView to prevent wrapping etc...
    const view = useMemo(() => {
        return config.dimension === "2d" ? new OrthographicView({
            id: `scatterplot-view-${id}`,
            controller: true,
            width,
            height,
            x: 0,
            y: 0,
        }) : new OrbitView({
            id: `scatterplot-view-${id}`,
            controller: true,
            width,
            height,
            x: 0,
            y: 0,
        });
    }, [width, height, config.dimension, id]);

    return (
        <>
        <DeckGL 
            layers={[scatterplotLayer]}
            useDevicePixels={true}
            // controller={true}
            viewState={viewState}
            // initialViewState={viewState} //consider not using react state for this
            views={view}
            onViewStateChange={v => {action(()=>config.viewState = v.viewState)()}}
        />
        </>
    );
});