import { ColorPaletteExtension, DetailView } from "@hms-dbmi/viv";
import { useCallback, useMemo, useRef } from "react";
import { shallow } from "zustand/shallow";
import { getFieldColor } from "../fieldColorManager";
import useGateLayers from "../hooks/useGateLayers";
import { useChartArrayGrid } from "../hooks/useChartArrayGrid";
import {
    useDensityGridCells,
    useDensityGridContours,
} from "../hooks/useDensityGridCells";
import { useVivDensityGridViewState } from "../hooks/useVivDensityGridViewState";
import { useChartID, useFilteredIndices } from "../hooks";
import { useSpatialLayers } from "../spatial_context";
import ChartArrayLayout from "./ChartArrayLayout";
import { cloneDeckLayer, getVivGridDetailViewId } from "./densityGridUtils";
import { useLoader, type OME_TIFF, useChannelsStore, useViewerStore, useViewerStoreApi } from "./avivatorish/state";
import MDVivViewer, { getVivId } from "./avivatorish/MDVivViewer";
import VivContrastExtension from "@/webgl/VivContrastExtension";
import { useOuterContainer } from "../screen_state";
import type { DeckGLProps, OrbitViewState, OrthographicViewState, PickingInfo } from "deck.gl";
import { getCombinedScatterTooltip } from "@/lib/scatterTooltip";
import { useOuterContainerDeckTooltip } from "../hooks/useOuterContainerDeckTooltip";

export default function VivDensityGridComponent() {
    const chartId = useChartID();
    const ome = useLoader() as OME_TIFF["data"];
    const viewerStore = useViewerStoreApi();
    const viewState = useViewerStore((store) => store.viewState);
    const outerContainer = useOuterContainer();
    const deckContainerRef = useRef<HTMLDivElement | null>(null);
    const rows = useFilteredIndices();

    const {
        scatterProps: { scatterplotLayer, greyScatterplotLayer, getTooltip, setScatterKeyboardActive },
        selectionLayer,
    } = useSpatialLayers();
    const { gateLabelLayer, gateDisplayLayer, controllerOptions } = useGateLayers();
    const { cells, densityFields, configuredFieldCount } = useDensityGridCells();

    const grid = useChartArrayGrid({
        chartId,
        cells,
        getViewId: (id, cell, index) => getVivGridDetailViewId(id, cell.key, index),
    });

    const contourLayers = useDensityGridContours(densityFields, grid.visibleCellIndices);
    const referenceCellWidth = grid.metrics.cellBounds[0]?.width ?? 0;
    const referenceCellHeight = grid.metrics.cellBounds[0]?.height ?? 0;
    const radiusPixels = useVivDensityGridViewState(referenceCellWidth, referenceCellHeight);

    const { colors, contrastLimits, channelsVisible, selections, brightness, contrast } = useChannelsStore(
        ({ colors, contrastLimits, channelsVisible, selections, brightness, contrast }) => ({
            colors,
            contrastLimits,
            channelsVisible,
            selections,
            brightness,
            contrast,
        }),
        shallow,
    );

    const extensions = useMemo(() => [new ColorPaletteExtension(), new VivContrastExtension()], []);
    const layerConfig = useMemo(
        () => ({
            loader: ome,
            selections,
            contrastLimits,
            extensions,
            colors,
            channelsVisible,
            brightness,
            contrast,
        }),
        [ome, selections, contrastLimits, extensions, colors, channelsVisible, brightness, contrast],
    );

    const detailViews = useMemo(
        () =>
            grid.visibleCellIndices
                .map((index) => {
                    const bounds = grid.getCellBounds(index);
                    const detailId = grid.viewIds[index];
                    if (!bounds || bounds.width <= 0 || bounds.height <= 0 || !detailId) return null;
                    return new DetailView({
                        id: detailId,
                        x: bounds.x,
                        y: bounds.y,
                        width: bounds.width,
                        height: bounds.height,
                        // @ts-expect-error viv runtime supports this, types do not
                        snapScaleBar: false,
                    });
                })
                .filter((view): view is DetailView => view !== null),
        [grid.visibleCellIndices, grid.viewIds, grid.metrics.cellBounds],
    );

    const viewStates = useMemo(() => {
        if (!viewState) return [];
        return grid.visibleCellIndices.map((index) => {
            const detailId = grid.viewIds[index];
            if (!detailId) return null;
            return { ...viewState, id: detailId };
        }).filter((vs): vs is NonNullable<typeof vs> => vs !== null);
    }, [grid.visibleCellIndices, grid.viewIds, viewState]);

    const deckLayers = useMemo(
        () =>
            grid.visibleCellIndices.flatMap((index) => {
                const contourLayer = contourLayers[index];
                const detailId = grid.viewIds[index];
                if (!contourLayer || !detailId) return [];

                const greyLayer = cloneDeckLayer(greyScatterplotLayer, {
                    id: `scatter-grey_${getVivId(detailId)}`,
                    viewId: detailId,
                });
                const scatterLayer = cloneDeckLayer(scatterplotLayer, {
                    id: `scatter_${getVivId(detailId)}`,
                    viewId: detailId,
                    contourLayers: [
                        {
                            ...contourLayer,
                            radiusPixels,
                            debounce: 250,
                            weightsTextureSize: 256,
                        },
                    ],
                });
                return [greyLayer, scatterLayer];
            }),
        [
            grid.visibleCellIndices,
            grid.viewIds,
            contourLayers,
            radiusPixels,
            scatterplotLayer,
            greyScatterplotLayer,
        ],
    );

    const gateLayers = useMemo(
        () => [gateDisplayLayer, gateLabelLayer].filter((layer) => layer !== null),
        [gateDisplayLayer, gateLabelLayer],
    );

    const allDeckLayers = useMemo(
        () =>
            selectionLayer
                ? [...deckLayers, ...gateLayers, selectionLayer]
                : [...deckLayers, ...gateLayers],
        [deckLayers, gateLayers, selectionLayer],
    );

    const getTooltipContent = useCallback(
        (info: PickingInfo) =>
            getCombinedScatterTooltip(info, {
                gateDisplayLayerId: gateDisplayLayer?.id,
                gateLabelLayerId: gateLabelLayer?.id,
                getPointTooltip: getTooltip,
            }),
        [gateDisplayLayer?.id, gateLabelLayer?.id, getTooltip],
    );

    const {
        clearTooltip,
        getTooltip: getPortalTooltip,
        suppressTooltipUntilPointerUp,
        tooltipPortal,
    } = useOuterContainerDeckTooltip(getTooltipContent, deckContainerRef);

    const deckProps: Partial<DeckGLProps> = useMemo(
        () => ({
            getTooltip: getPortalTooltip,
            layers: allDeckLayers,
            controller: {
                doubleClickZoom: false,
                dragPan: controllerOptions.dragPan,
            },
        }),
        [allDeckLayers, getPortalTooltip, controllerOptions.dragPan],
    );

    const renderCell = useCallback(
        (index: number) => {
            const field = densityFields[index];
            if (!field) return null;
            const color = getFieldColor(field.field);
            return (
                <div
                    className="pointer-events-none absolute left-2 top-2 z-10 max-w-[calc(100%-16px)] truncate rounded-sm px-2 py-1 text-xs"
                    style={{
                        backgroundColor: "rgba(15, 23, 42, 0.76)",
                        color: `rgb(${color[0]}, ${color[1]}, ${color[2]})`,
                    }}
                    title={field.name || field.field}
                >
                    {field.name || field.field}
                </div>
            );
        },
        [densityFields],
    );

    const hasCanvas =
        grid.hasCanvas &&
        !!ome &&
        !!viewState &&
        Number.isFinite(radiusPixels) &&
        radiusPixels > 0 &&
        detailViews.length > 0;

    const vivOverlay = useMemo(
        () =>
            hasCanvas ? (
                <MDVivViewer
                    outerContainer={outerContainer}
                    selectionLayer={selectionLayer}
                    views={detailViews}
                    layerProps={detailViews.map(() => layerConfig)}
                    viewStates={viewStates}
                    useDevicePixels={true}
                    onViewStateChange={(e: {
                        viewId: string;
                        viewState: OrthographicViewState | OrbitViewState;
                    }) => {
                        viewerStore.setState({
                            viewState: { ...e.viewState, id: e.viewId },
                        });
                    }}
                    deckProps={deckProps}
                />
            ) : null,
        [hasCanvas, outerContainer, selectionLayer, detailViews, layerConfig, viewStates, viewerStore, deckProps],
    );

    if (!viewState) {
        return <div className="flex h-full items-center justify-center text-sm">Loading...</div>;
    }
    if (configuredFieldCount === 0) {
        return <div className="flex h-full items-center justify-center text-sm">Choose density fields to build the grid.</div>;
    }
    if (densityFields.length === 0) {
        return <div className="flex h-full items-center justify-center text-sm">Loading density fields...</div>;
    }
    if (rows.length === 0) {
        return <div className="flex h-full items-center justify-center text-sm">No rows remain after the current filters.</div>;
    }

    return (
        <>
            <div
                ref={deckContainerRef}
                className="relative flex h-full w-full min-w-0 flex-col"
                aria-label="Spatial scatter density grid"
                style={{ outline: "none" }}
                onPointerDown={suppressTooltipUntilPointerUp}
                onMouseDown={() => setScatterKeyboardActive(true)}
                onMouseEnter={() => setScatterKeyboardActive(true)}
                onMouseLeave={() => {
                    clearTooltip();
                    setScatterKeyboardActive(false);
                }}
            >
                <div
                    ref={grid.scrollContainerRef}
                    className="min-h-0 min-w-0 flex-1 overflow-auto overflow-x-hidden"
                >
                    <ChartArrayLayout
                        cellCount={grid.cellCount}
                        cellKeys={grid.cellKeys}
                        layoutRef={grid.layoutRef}
                        canvasOverlay={vivOverlay}
                        renderCell={renderCell}
                    />
                </div>
            </div>
            {tooltipPortal}
        </>
    );
}
