import DeckGL from "@deck.gl/react";
import { OrthographicView, OrbitView } from '@deck.gl/core';
import { observer } from "mobx-react-lite";
import { useChartSize, useConfig, useFilterArray, useFilteredIndices, useParamColumns, useSimplerFilteredIndices } from "../hooks";
import { ScatterplotLayer } from "@deck.gl/layers";
import { DataFilterExtension } from "@deck.gl/extensions";
import { useEffect, useId, useMemo, useState } from "react";
import { useChart, useDataStore } from "../context";
import type { DeckScatterConfig } from "./DeckScatterReactWrapper";
import { action } from "mobx";
import type { OrbitViewState } from "@deck.gl/core";
import type { LoadedDataColumn } from "@/charts/charts";
import { Axis, Scale } from "@visx/visx";
import "../../charts/css/charts.css";
import { SpatialAnnotationProvider, useSpatialLayers } from "../spatial_context";
import SelectionOverlay from "./SelectionOverlay";

const margin = { top: 10, right: 10, bottom: 40, left: 60 };
/** todo this should be common for viv / scatter_state, pending refactor
 * there should be hooks getting the range of a filtered column 
 * so that multiple charts can re-use the computation (but if not used, it's not computed)
 */
function useZoomOnFilter(data: Uint32Array) {
    const [width, height] = useChartSize();
    const [cx, cy, cz] = useParamColumns();
    // const data = useFilteredIndices(); //<< does calling this multiple times mean wasting space?
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
        // ^ maybe we could have properties on the columns that are the min/max of the filtered data
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
                zoom: Math.log2(Math.min(width / (maxX - minX), height / (maxY - minY))) - 0.2,
                // any kind of interpolator would probably mean changing how we handle viewState
                // transitionInterpolator: new FlyToInterpolator({speed: 1.2}),
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


/**
 * Currently this implementation is somewhat separate from VivMDVReact / scatter_state,
 * but it should be more unified in the future.
 * 
 * Things that need work:
 * - GeoJsonEditableLayer & associated selection state
 *  - think about selection in 3d
 * - density (in future, more different types of layers)
 *  - fancy ways of e.g. selecting things based on density (metaballs)
 * - fix mouse events in popout
 * - auto-scaling of radius (including when parameter changes)
 * - different scaling along different axes when not using physical units
 * - highlighting
 * - tooltips
 *   - enhanced version with more info (in general for all charts, or at least new ones)
 * - axis configuration (including via direct manipulation)
 */
const DeckScatter = observer(function DeckScatterComponent() {
    const id = useId();
    const [width, height] = useChartSize();
    const [cx, cy, cz] = useParamColumns() as LoadedDataColumn<"double">[];
    // const data = useSimplerFilteredIndices();
    const data = useFilteredIndices(); //changed to fallback to simplerFilteredIndices when filterColumn is not set
    const config = useConfig<DeckScatterConfig>();
    const { opacity, radius, course_radius, viewState, on_filter, color_by } = config;
    // const is2d = dimension === "2d";
    //todo more clarity on radius units - but large radius was causing big problems after deck upgrade
    const radiusScale = radius * course_radius * 0.01;// * (is2d ? 1 : 0.01);
    const chart = useChart();
    //todo colorBy should be done differently (also bearing in mind multiple layers)
    const colorBy = (chart as any).colorBy;
    // const colorBy = color_by;
    const chartWidth = width - margin.left - margin.right;
    const chartHeight = height - margin.top - margin.bottom;
    useZoomOnFilter(data);

    const greyOnFilter = on_filter === "grey";

    const { scatterProps, selectionLayer } = useSpatialLayers();
    // this is now somewhat able to render for "2d", pending further tweaks
    const { scatterplotLayer, getTooltip } = scatterProps;


    const XscatterplotLayer = useMemo(() => new ScatterplotLayer({
        id: `scatterplot-layer-${id}`,
        data,
        pickable: true,
        opacity,
        stroked: false,
        filled: true,
        radiusScale,
        getPosition: (index, { target }) => {
            target[0] = cx.data[index];
            target[1] = cy.data[index];
            if (cz) target[2] = cz.data[index];
            return target as [number, number];
        },
        // maybe we want a kind of concrete config synthesis that means we have colorBy as a function...
        // also maybe getPosition...
        // to what extent is that stuff we want to generalise?
        // ! note if we want to have `data` include all indices & filter with `DataFilterExtension`,
        // (which would help with transitions on filtered data)
        // it ends up being much slower than just using the filtered indices
        // and we need to adapt the colorBy function as it will no longer be passed the index as first argument.
        getFillColor: colorBy ?? [61, 126, 180],
        getLineColor: [0, 0, 0],
        updateTriggers: {
            getFillColor: colorBy,
            getPosition: [cx.data, cy.data, cz?.data],
        },
        billboard: true,
        parameters: {
            depthTest: false,
        },
        transitions: {
            //todo make this common for all layers, control from config...
            getPosition: {
                duration: 100,
                //https://easings.net/#easeInOutCubic
                easing: (x: number) => x < 0.5 ? 4 * x * x * x : 1 - (-2 * x + 2) ** 3 / 2,
                // type: "spring",
            },
        }
    }), [data, cx, cy, cz, colorBy, opacity, radiusScale, id]);

    const filterValue = useFilterArray();

    const greyScatterplotLayer = useMemo(() => new ScatterplotLayer({
        id: `scatterplot-layer-grey-${id}`,
        data: { length: cx.data.length },
        // pickable: true,
        opacity,
        stroked: false,
        filled: true,
        radiusScale,
        getPosition: (_, { target, index }) => {
            target[0] = cx.data[index];
            target[1] = cy.data[index];
            if (cz) target[2] = cz.data[index];
            return target as [number, number];
        },
        getFillColor: [200, 200, 200],
        getLineColor: [0, 0, 0],
        billboard: true,
        parameters: {
            depthTest: false,
        },
        transitions: {
            getPosition: {
                duration: 100,
                //https://easings.net/#easeInOutCubic
                easing: (x: number) => x < 0.5 ? 4 * x * x * x : 1 - (-2 * x + 2) ** 3 / 2,
                // type: "spring",
            },
        },
        getFilterValue: (_: any, { index }: { index: number }) => filterValue[index] || 0, //how do we just pass the buffer?
        filterRange: [0.5, 1],
        updateTriggers: {
            //! using `data` as a trigger as `filterValue` was behind by one frame or something
            getFilterValue: data,
            getPosition: [cx.data, cy.data, cz?.data],
        },
        extensions: [new DataFilterExtension()],
        visible: greyOnFilter,
    }), [cx, cy, cz, opacity, radiusScale, id, filterValue, data, greyOnFilter]);
    
    // we need an OrthographicView to prevent wrapping etc...
    const view = useMemo(() => {
        return config.dimension === "2d" ? new OrthographicView({
            id: `scatterplot-view-${id}`,
            controller: true,
            width: chartWidth,
            height: chartHeight,
            x: margin.left + 1,
            y: margin.top - 1,
            flipY: false,
        }) : new OrbitView({
            id: `scatterplot-view-${id}`,
            controller: true,
            width: chartWidth,
            height: chartHeight,
            x: margin.left + 1,
            y: margin.top - 1,
        });
    }, [chartWidth, chartHeight, config.dimension, id]);
    //! deck doesn't like it if we change the layers array - better to toggle visibility
    const layers = [selectionLayer, scatterplotLayer, greyScatterplotLayer];
    // todo - these need to be encapsulated better, the DeckGL component should be in a smaller
    // area with the axes outside of it.
    // axes need to respond to the viewState...
    // This implementation will mean that the whole react component will re-render when the viewState changes.
    // This is not very good for performance - we may consider using refs or something to avoid this, 
    // and/or debouncing the viewState changes.
    const ranges = useMemo(() => {
        viewState;
        // first time around, we get an exception because scatterplotLayer hasn't been rendered yet
        try {
            const p = scatterplotLayer.unproject([0, 0]);
            const p2 = scatterplotLayer.unproject([chartWidth, chartHeight]);
            const domainX = [p[0], p2[0]];
            const domainY = [p2[1], p[1]];
            return { domainX, domainY };
        } catch (e) {
            return { domainX: cx.minMax, domainY: cy.minMax };
        }
    }, [cx.minMax, cy.minMax, chartWidth, chartHeight, scatterplotLayer, viewState]);
    const scaleX = useMemo(() => Scale.scaleLinear({
        domain: ranges.domainX, // e.g. [min, max]
        range: [margin.left, chartWidth + margin.left],
    }), [chartWidth, ranges]);
    const scaleY = useMemo(() => Scale.scaleLinear({
        domain: ranges.domainY, // e.g. [min, max]
        range: [chartHeight + margin.top, margin.top],
    }), [chartHeight, ranges]);

    return (
        <>
            <DeckGL
                layers={layers}
                useDevicePixels={true}
                controller={true}
                viewState={viewState}
                // initialViewState={viewState} //consider not using react state for this        
                views={view}
                onViewStateChange={v => { action(() => config.viewState = v.viewState)() }}
            />
            <svg width={width} height={height}>
                <Axis.AxisBottom
                    top={chartHeight + margin.top}
                    scale={scaleX}
                    stroke={"var(--text_color)"}
                    tickStroke={"var(--text_color)"}
                    tickLabelProps={() => ({
                        fill: "var(--text_color)",
                        fontSize: 9.5,
                        textAnchor: "middle",
                    })}
                    labelProps={{
                        fill: "var(--text_color)",
                        fontSize: 10,
                        textAnchor: "middle",
                    }}
                    label={cx.name}
                />
                <Axis.AxisLeft
                    left={margin.left}
                    scale={scaleY}
                    stroke={"var(--text_color)"}
                    tickStroke={"var(--text_color)"}
                    tickLabelProps={() => ({
                        fill: "var(--text_color)",
                        fontSize: 9.5,
                        textAnchor: "end",
                    })}
                    labelProps={{
                        fill: "var(--text_color)",
                        fontSize: 10,
                    }}
                    label={cy.name}
                />
            </svg>
        </>
    );
});

export default () => {
    const chart = useChart();
    // in order for SelectionOverlay to work, we need to review how our implementation works
    // vs useScatterplotLayer in spatial_context.tsx
    return (
        <SpatialAnnotationProvider chart={chart}>
            <SelectionOverlay />
            <DeckScatter />
        </SpatialAnnotationProvider>
    );
};