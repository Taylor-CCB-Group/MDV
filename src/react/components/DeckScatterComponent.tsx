import DeckGL from "@deck.gl/react";
import { OrthographicView, OrbitView } from '@deck.gl/core';
import { observer } from "mobx-react-lite";
import { useChartSize, useConfig, useFilterArray, useFilteredIndices, useParamColumns } from "../hooks";
import { GeoJsonLayer, LineLayer, ScatterplotLayer, TextLayer } from "@deck.gl/layers";
import { DataFilterExtension } from "@deck.gl/extensions";
import { useCallback, useEffect, useId, useMemo, useRef } from "react";
import { useChart } from "../context";
import type { DeckScatterConfig } from "./DeckScatterReactWrapper";
import { action } from "mobx";
import type { DataColumn, LoadedDataColumn, NumberDataType } from "@/charts/charts";
import { allNumeric } from "@/lib/columnTypeHelpers";
import { SpatialAnnotationProvider, useSpatialLayers } from "../spatial_context";
import SelectionOverlay from "./SelectionOverlay";
import { useScatterRadius } from "../scatter_state";
import AxisComponent from "./AxisComponent";
import { useOuterContainer } from "../screen_state";
import { rebindMouseEvents } from "@/lib/deckMonkeypatch";
import { useGateManager } from "../gates/useGateManager";
import { computeCentroid } from "../gates/gateUtils";

//todo this should be in a common place etc.
const colMid = ({minMax}: DataColumn<NumberDataType>) => minMax[0] + (minMax[1] - minMax[0]) / 2;

// const margin = { top: 10, right: 10, bottom: 40, left: 60 };
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
        action(() => {
            chart.ignoreStateUpdate = true;
            chart.pendingRecenter = false;
            
            // Calculate range
            const dx = maxX - minX;
            const dy = maxY - minY;

            // Use relative epsilon based on data scale or absolute minimum
            const epsilon = Math.max(1e-9, Math.min(Math.abs(minX), Math.abs(maxX)) * 1e-6);
            const safeDx = Math.max(epsilon, dx);
            const safeDy = Math.max(epsilon, dy);

            // Calculate zoom with fallback and NaN check
            let zoom = Math.log2(Math.min(width / safeDx, height / safeDy)) - 0.6;
            if (!Number.isFinite(zoom)) {
                console.warn("Invalid zoom value detected, using fallback");
                zoom = 0; // Default zoom fallback
            }

            config.viewState = {
                target: [(minX + maxX) / 2, (minY + maxY) / 2, 0],
                zoom,
                // any kind of interpolator would probably mean changing how we handle viewState
                // transitionInterpolator: new FlyToInterpolator({speed: 1.2}),
            };
            if (config.dimension === "3d") {
                const vs = config.viewState;
                vs.rotationOrbit = vs.rotationOrbit || 0;
                vs.rotationX = vs.rotationX || 0;
                if (!config.zoom_on_filter) {
                    if (!cz) {
                        console.warn("3D viewState requires z coordinate");
                        return;
                    }
                    const params = [cx, cy, cz];
                    if (!allNumeric(params)) {
                        throw new Error("3D viewState requires numeric coordinates");
                    }
                    vs.target = params.map(colMid) as [number, number, number];
                }
            }

            chart.ignoreStateUpdate = false;
        })();
    }, [data, cx, cy, cz, width, height, config, pendingRecenter, chart]);
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
    const [cx, cy, ...density] = useParamColumns() as LoadedDataColumn<"double">[];
    const cz = useParamColumns()[2] as LoadedDataColumn<"double">;
    const data = useFilteredIndices(); //changed to fallback to simplerFilteredIndices when filterColumn is not set
    const config = useConfig<DeckScatterConfig>();
    const { opacity, viewState, on_filter, dimension } = config;
    const is2d = dimension === "2d";
    //todo more clarity on radius units - but large radius was causing big problems after deck upgrade
    const radiusScale = useScatterRadius();
    //todo colorBy should be done differently (also bearing in mind multiple layers)
    
    //todo this shouldn't be repeated here and in AxisComponent
    const xSize = is2d ? config.axis.x.size : 0;
    const ySize = is2d ? config.axis.y.size : 0;

    const margin = useMemo(() => (
        is2d ? {
            top: 10,
            right: 10,
            bottom: xSize,
            left: ySize,
        } : {
            top: 0,
            right: 0,
            bottom: 0,
            left: 0,
        }
    ), [is2d, xSize, ySize]);
    const chartWidth = width - margin.left - margin.right;
    //there could be a potential off-by-one/two error somewhere down the line
    //if we don't fully understand reasons for `- 3.5` here.
    //prevents overlapping with x-axis.
    const chartHeight = height - margin.top - margin.bottom - 3.5;
    useZoomOnFilter(data);

    const greyOnFilter = on_filter === "grey";

    const { scatterProps, selectionLayer } = useSpatialLayers();
    // this is now somewhat able to render for "2d", pending further tweaks
    //! beware unproject from here is not what we want, should review
    const { scatterplotLayer, getTooltip } = scatterProps;
    const gateManager = useGateManager();

    const filterValue = useFilterArray();

    // this should move in to scatter_state, common with viv...
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

    const axisLinesLayer = useMemo(() => {
        if (is2d) return null;
        return new LineLayer({
            id: `axis-lines-${id}`,
            data: [
                { sourcePosition: [0, 0, 0], targetPosition: [cx.minMax[1], 0, 0], color: [255, 0, 0] },
                { sourcePosition: [0, 0, 0], targetPosition: [0, cy.minMax[1], 0], color: [0, 255, 0] },
                { sourcePosition: [0, 0, 0], targetPosition: [0, 0, cz.minMax[1]], color: [0, 0, 255] },
            ],
            getSourcePosition: (d: any) => d.sourcePosition,
            getTargetPosition: (d: any) => d.targetPosition,
            getColor: (d: any) => d.color,
            getWidth: 1,
            updateTriggers: {
                getSourcePosition: [cx.minMax, cy.minMax, cz.minMax],
                getTargetPosition: [cx.minMax, cy.minMax, cz.minMax],
            },
        });
    }, [is2d, id, cx.minMax, cy.minMax, cz?.minMax]);
    
    // we need an OrthographicView to prevent wrapping etc...
    // if in future we have subgraphs sharing a canvas, we will need to
    // make sure that the view is set up correctly for each subgraph.
    const view = useMemo(() => {
        return config.dimension === "2d" ? new OrthographicView({
            id: `scatterplot-view-${id}`,
            controller: true,
            width: chartWidth,
            height: chartHeight,
            x: 0,
            y: 0,
            flipY: false,
        }) : new OrbitView({
            id: `scatterplot-view-${id}`,
            controller: true,
            width: chartWidth,
            height: chartHeight,
            x: 0,
            y: 0,
        });
    }, [chartWidth, chartHeight, config.dimension, id]);

    // Gate labels layer - access gatesArray to ensure MobX reactivity
    const allGates = gateManager.gatesArray;
    const gateLabelsLayer = useMemo(() => {
        if (!cx || !cy || allGates.length === 0) return null;
        
        // Get gates that match the current chart's X/Y columns
        const relevantGates = allGates.filter(
            gate => gate.columns[0] === cx.field && gate.columns[1] === cy.field
        );
        
        if (relevantGates.length === 0) return null;
        
        // Create data array with gate positions and names
        const gateData = relevantGates.map(gate => {
            const centroid = computeCentroid(gate.geometry);
            // For 2D, use z=0; for 3D, use middle of z range if cz exists
            const z = is2d ? 0 : (cz?.minMax ? (cz.minMax[0] + cz.minMax[1]) / 2 : 0);
            return {
                position: [centroid[0], centroid[1], z] as [number, number, number],
                name: gate.name,
            };
        });
        
        return new TextLayer({
            id: `gate-labels-${id}`,
            data: gateData,
            getPosition: (d: { position: [number, number, number] }) => d.position,
            getText: (d: { name: string }) => d.name,
            getSize: 12,
            getColor: [255, 255, 255, 255],
            getTextAnchor: 'middle',
            getAlignmentBaseline: 'center',
            getBackgroundColor: [0, 0, 0, 191], // Semi-transparent black background
            getBorderColor: [255, 255, 255, 76], // Semi-transparent white border
            getBorderWidth: 1,
            padding: [3, 8],
            billboard: true,
            pickable: false,
            updateTriggers: {
                getPosition: [allGates.length, relevantGates.map(g => g.id).join(',')],
                getText: [allGates.length, relevantGates.map(g => g.name).join(',')],
            },
        });
    }, [id, is2d, cx, cy, cz, allGates]);
    
    // Gate polygons layer - renders the actual gate shapes
    const gatePolygonsLayer = useMemo(() => {
        if (!cx || !cy || allGates.length === 0) return null;
        
        // Get gates that match the current chart's X/Y columns
        const relevantGates = allGates.filter(
            gate => gate.columns[0] === cx.field && gate.columns[1] === cy.field
        );
        
        if (relevantGates.length === 0) return null;
        
        // Combine all gate geometries into a single FeatureCollection
        // Each gate's features are copied with gate metadata in properties
        const combinedFeatures = relevantGates.flatMap(gate => 
            gate.geometry.features.map(feature => ({
                ...feature,
                properties: {
                    ...feature.properties,
                    gateId: gate.id,
                    gateName: gate.name
                }
            }))
        );
        
        return new GeoJsonLayer({
            id: `gate-polygons-${id}`,
            data: {
                type: "FeatureCollection",
                features: combinedFeatures
            } as any,
            filled: true,
            getFillColor: [76, 175, 80, 30], // Semi-transparent green fill
            getLineColor: [76, 175, 80, 200], // Green border (more opaque)
            getLineWidth: 2,
            lineWidthMinPixels: 1,
            pickable: false, // Gates are not interactive in MVP
            updateTriggers: {
                updateTriggers: {
                    getFillColor: [allGates.length, relevantGates.map(g => g.id).join(',')],
                    getLineColor: [allGates.length, relevantGates.map(g => g.id).join(',')],
                },
            },
        });
    }, [id, cx, cy, allGates]);
    //! deck doesn't like it if we change the layers array - better to toggle visibility
    const layers = [scatterplotLayer, greyScatterplotLayer, gatePolygonsLayer, selectionLayer, axisLinesLayer, 
         gateLabelsLayer
    ].filter(x => x !== null);
    
    // unproject used for updating ranges - may refactor hooks around this
    const unproject = useCallback((coords: [number, number]) => {
        // make sure it applies to the right `this`
        return scatterplotLayer.unproject(coords);
    }, [scatterplotLayer]);
    const outerContainer = useOuterContainer();
    const deckRef = useRef<any>();

    useEffect(() => {
        // useParamColumns() only returns when columns are loaded (it throws otherwise)
        // So if we get here, columns are ready
        if (cx && cy) {
            gateManager.updateGatesColumnWhenReady().catch(console.error);
        }
    }, [cx, cy, gateManager]);

    // biome-ignore lint/correctness/useExhaustiveDependencies: selectionLayer might change without us caring
    useEffect(() => {
        outerContainer; // make sure the hook runs when this changes
        if (deckRef.current) {
            try {
                // const deck: Deck<any> = deckRef.current.deck;
                const deck = deckRef.current.deck;// as Deck<any>;
                return rebindMouseEvents(deck, selectionLayer);
            } catch (e) {
                console.error(
                    "attempt to reset deck eventManager element failed - could be related to brittle deck monkeypatch",
                    e,
                );
            }
        }
    }, [outerContainer]);
    
    // we want default controller options, but we want a new one when the outerContainer changes
    // this doesn't seem to help re-register mouse events.
    // const controller = useMemo(() => ({inertia: 10+Math.random()}), [outerContainer])
    return (
        <>
            <AxisComponent config={config} unproject={unproject}>
                <DeckGL
                    ref={deckRef}
                    layers={layers}
                    useDevicePixels={true}
                    controller={true}
                    viewState={viewState}
                    // initialViewState={viewState} //consider not using react state for this        
                    views={view}
                    onViewStateChange={v => { action(() => config.viewState = v.viewState)() }}
                    getTooltip={getTooltip}
                />
            </AxisComponent>
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