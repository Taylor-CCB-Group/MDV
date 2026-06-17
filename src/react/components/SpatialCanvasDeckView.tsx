import {
    SpatialViewer,
    useSpatialCanvasRenderer,
    type ViewState as SpatialCanvasViewState,
} from "@spatialdata/vis";
import type { SpatialData } from "@spatialdata/core";
import type { Layer } from "@deck.gl/core";
import type { DeckGLProps } from "deck.gl";
import { observer } from "mobx-react-lite";
import { useEffect, useMemo } from "react";

import { useChart, useDataStore } from "../context";
import { useConfig } from "../hooks";
import { inferTableAssociation } from "@/react/spatial_table_association";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "./SpatialDataMDVReact";
import { useViewerStoreApi } from "./avivatorish/state";
import type { OrbitViewState, OrthographicViewState } from "deck.gl";

export const SpatialCanvasDeckView = observer(function SpatialCanvasDeckView({
    spatialData,
    coordinateSystem,
    spatialdataPath,
    spatialViewState,
    width,
    height,
    deckLayers,
    unifiedLayerOrder,
    deckProps,
    toMdvViewState,
}: {
    spatialData: SpatialData | null;
    coordinateSystem: string | null;
    spatialdataPath: string | null;
    spatialViewState: SpatialCanvasViewState | null;
    width: number;
    height: number;
    deckLayers: Layer[];
    unifiedLayerOrder: string[];
    deckProps: Partial<DeckGLProps>;
    toMdvViewState: (
        next: SpatialCanvasViewState,
        current: OrthographicViewState | OrbitViewState | null | undefined,
    ) => OrthographicViewState | OrbitViewState;
}) {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const config = useConfig<SpatialDataMdvReactConfig>();
    const dataStore = useDataStore();
    const viewerStore = useViewerStoreApi();

    // Subscribe to lightweight canvas patches (opacity, contrast, etc.).
    void chart.canvasLayerRevision;

    const { layers, layerOrder } = chart.canvasLayerState;

    const renderer = useSpatialCanvasRenderer({
        spatialData,
        coordinateSystem,
        layers,
        layerOrder,
        viewState: spatialViewState,
        onViewStateChange: (nextViewState) => {
            viewerStore.setState({
                viewState: toMdvViewState(
                    nextViewState,
                    viewerStore.getState().viewState,
                ),
            });
        },
        width,
        height,
        deckLayers,
        autoFit: true,
    });

    const rendererContextValue = useMemo(
        () => ({
            spatialData,
            coordinateSystem,
            availableElements: renderer.availableElements,
            getImageLayerLoadedData: renderer.getImageLayerLoadedData,
            getLabelsLayerLoadedData: renderer.getLabelsLayerLoadedData,
            getLayerLoadState: renderer.getLayerLoadState,
            isLoading: renderer.isLoading,
            isBlocking: renderer.isBlocking,
            inferTableAssociation: (elementKey?: string | null) =>
                inferTableAssociation({
                    spatialData,
                    spatialdataPath,
                    regionId: config.region,
                    elementKey,
                    dataStore,
                }),
        }),
        [
            coordinateSystem,
            config.region,
            dataStore,
            renderer.availableElements,
            renderer.getImageLayerLoadedData,
            renderer.getLabelsLayerLoadedData,
            renderer.getLayerLoadState,
            renderer.isBlocking,
            renderer.isLoading,
            spatialData,
            spatialdataPath,
        ],
    );

    useEffect(() => {
        chart.spatialRendererContext = rendererContextValue;
    }, [chart, rendererContextValue]);

    return (
        <>
            {renderer.isBlocking && (
                <div className="absolute right-2 top-2 rounded bg-black/70 px-2 py-1 text-[11px] text-white">
                    Loading layer data...
                </div>
            )}
            {!renderer.isBlocking && renderer.isLoading && (
                <div className="absolute right-2 top-2 rounded bg-neutral-900/80 px-2 py-1 text-[11px] text-white">
                    Refreshing layer metadata...
                </div>
            )}
            <SpatialViewer
                width={width}
                height={height}
                viewState={spatialViewState}
                onViewStateChange={(nextViewState) => {
                    viewerStore.setState({
                        viewState: toMdvViewState(
                            nextViewState,
                            viewerStore.getState().viewState,
                        ),
                    });
                }}
                layers={renderer.deckLayers}
                layerOrder={unifiedLayerOrder}
                vivLayerProps={
                    renderer.vivLayerProps.length > 0
                        ? renderer.vivLayerProps
                        : undefined
                }
                deckProps={deckProps}
            />
        </>
    );
});
