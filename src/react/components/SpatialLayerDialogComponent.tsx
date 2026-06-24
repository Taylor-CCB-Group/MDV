import ArrowDropDownIcon from "@mui/icons-material/ArrowDropDown";
import DragIndicatorIcon from "@mui/icons-material/DragIndicator";
import HighlightOffIcon from "@mui/icons-material/HighlightOff";
import {
    Accordion,
    AccordionDetails,
    AccordionSummary,
    Autocomplete,
    Checkbox,
    IconButton,
    Slider,
    TextField,
    Typography,
} from "@mui/material";
import {
    DndContext,
    PointerSensor,
    closestCenter,
    useSensor,
    useSensors,
    type DragEndEvent,
} from "@dnd-kit/core";
import {
    SortableContext,
    useSortable,
    verticalListSortingStrategy,
} from "@dnd-kit/sortable";
import { CSS } from "@dnd-kit/utilities";
import type { RenderStackSpatialElementType } from "@spatialdata/layers";
import type { LayerConfig } from "@spatialdata/vis";
import { useSpatialData } from "@spatialdata/react";
import { observer } from "mobx-react-lite";
import { useCallback, useEffect, useMemo } from "react";

import {
    DECK_OVERLAY_IDS,
    DECK_OVERLAY_LABELS,
    deckHostLayerId,
    deckIdFromHostLayerId,
    type DeckOverlayId,
} from "@/react/spatialdata/host_overlay_ids";
import { renderStackEntryDisplayName, renderStackOrderLabel } from "@/react/spatialdata/render_stack_display";
import {
    insertHostRenderStackEntry,
    insertSpatialRenderStackEntry,
    isRemovableRenderStackEntry,
    removeRenderStackEntry,
    reorderRenderStackEntries,
    renderStackEntryIds,
    seedRenderStackFromSpatialData,
    useRenderStackEntry,
    useRenderStackMutation,
} from "@/react/spatialdata/render_stack_control";
import {
    defaultPropsForSpatialElement,
    listAvailableSpatialEntries,
} from "@/react/spatialdata/render_stack_defaults";
import { NO_TABLE_ASSOCIATION } from "@/react/spatialdata/table_association";
import { useImageLayerPanelDefaults } from "@/react/spatialdata/use_image_layer_panel_defaults";
import { useChart, useDataStore } from "../context";
import { useConfig } from "../hooks";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "./SpatialDataMDVReact";
import DeckOverlayLayerPanel from "./spatialLayers/DeckOverlayLayerPanel";
import ImageLayerPanel from "./spatialLayers/ImageLayerPanel";
import LabelsLayerPanel from "./spatialLayers/LabelsLayerPanel";
import PointsLayerPanel from "./spatialLayers/PointsLayerPanel";
import ShapesLayerPanel from "./spatialLayers/ShapesLayerPanel";

type InsertOption =
    | { kind: "spatial"; type: RenderStackSpatialElementType; elementKey: string; label: string }
    | { kind: "host"; deckId: DeckOverlayId; label: string };

function getAvailableFields(dataStore: ReturnType<typeof useDataStore>): string[] {
    return Object.keys(dataStore.columnIndex).sort();
}

const LayerDetails = observer(function LayerDetails({ entryId }: { entryId: string }) {
    const dataStore = useDataStore();
    const chartConfig = useConfig<SpatialDataMdvReactConfig>();
    const { spatialData } = useSpatialData();
    const { entry, layer, patchLayer } = useRenderStackEntry(entryId);
    const availableFields = useMemo(() => getAvailableFields(dataStore), [dataStore]);
    const chartColorBy =
        typeof chartConfig.color_by === "string" ? chartConfig.color_by : undefined;
    const elementKey = entry?.kind === "spatial" ? entry.source.elementKey : undefined;
    const { channelNames, loaderDefaults } = useImageLayerPanelDefaults(spatialData, elementKey);

    if (!entry) return null;

    if (entry.kind === "host") {
        const deckId = deckIdFromHostLayerId(entry.source.hostLayerId);
        if (!deckId) return null;
        return <DeckOverlayLayerPanel deckId={deckId} />;
    }

    if (entry.kind !== "spatial" || !layer) return null;

    switch (entry.source.elementType) {
        case "image":
            return (
                <ImageLayerPanel
                    entryId={entryId}
                    channelNames={channelNames}
                    loaderDefaults={loaderDefaults}
                />
            );
        case "shapes":
            return (
                <ShapesLayerPanel
                    config={layer as Extract<LayerConfig, { type: "shapes" }>}
                    association={NO_TABLE_ASSOCIATION}
                    availableFields={availableFields}
                    chartColorBy={chartColorBy}
                    updateLayer={patchLayer}
                />
            );
        case "points":
            return (
                <PointsLayerPanel
                    config={layer as Extract<LayerConfig, { type: "points" }>}
                    updateLayer={patchLayer}
                />
            );
        case "labels":
            return (
                <LabelsLayerPanel
                    config={layer as Extract<LayerConfig, { type: "labels" }>}
                    association={NO_TABLE_ASSOCIATION}
                    availableFields={availableFields}
                    updateLayer={patchLayer}
                />
            );
        default:
            return null;
    }
});

const SortableLayerAccordion = observer(function SortableLayerAccordion({
    entryId,
    onRemove,
}: {
    entryId: string;
    onRemove: (entryId: string) => void;
}) {
    const { entry, patchEntry, patchProps } = useRenderStackEntry(entryId);
    const {
        attributes,
        listeners,
        setNodeRef,
        transform,
        transition,
        isDragging,
    } = useSortable({ id: entryId });

    if (!entry) return null;

    const opacity = typeof entry.props.opacity === "number" ? entry.props.opacity : 1;
    const supportsOpacity = entry.kind === "spatial";

    const style = {
        transform: CSS.Translate.toString(transform),
        transition,
        opacity: isDragging ? 0.5 : 1,
    };

    return (
        <Accordion
            ref={setNodeRef}
            style={style}
            disableGutters
            sx={{ width: "100%", display: "block" }}
            defaultExpanded={entry.kind === "spatial"}
        >
            <AccordionSummary expandIcon={<ArrowDropDownIcon />}>
                {isRemovableRenderStackEntry(entry) && (
                    <IconButton
                        onClick={(event) => {
                            event.stopPropagation();
                            onRemove(entryId);
                        }}
                        aria-label="remove layer"
                        size="small"
                        sx={{
                            position: "absolute",
                            right: "-18px",
                            top: "-18px",
                        }}
                    >
                        <HighlightOffIcon fontSize="small" />
                    </IconButton>
                )}
                <div className="flex w-full items-center gap-2">
                    <IconButton
                        {...attributes}
                        {...listeners}
                        size="small"
                        sx={{ cursor: "grab" }}
                    >
                        <DragIndicatorIcon fontSize="small" />
                    </IconButton>
                    <Checkbox
                        size="small"
                        checked={entry.visible}
                        onClick={(event) => event.stopPropagation()}
                        onChange={(event) =>
                            patchEntry({ visible: event.target.checked })
                        }
                    />
                    <Typography variant="subtitle1" className="grow">
                        {renderStackEntryDisplayName(entry)}
                    </Typography>
                    {supportsOpacity && (
                        <div
                            className="w-32 px-2"
                            onClick={(event) => event.stopPropagation()}
                        >
                            <Slider
                                size="small"
                                min={0}
                                max={1}
                                step={0.01}
                                value={opacity}
                                onChange={(_, value) => {
                                    if (typeof value === "number") {
                                        patchProps({ opacity: value });
                                    }
                                }}
                            />
                        </div>
                    )}
                </div>
            </AccordionSummary>
            <AccordionDetails>
                <LayerDetails entryId={entryId} />
            </AccordionDetails>
        </Accordion>
    );
});

const SpatialLayerDialogComponent = observer(() => {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const config = useConfig<SpatialDataMdvReactConfig>();
    const { spatialData } = useSpatialData();
    const stack = config.renderStack;
    const mutateStack = useRenderStackMutation();
    const coordinateSystem = //todo datasources schema including spatial regions when referring to dataStore here.
        (chart.dataStore.regions?.all_regions[config.region] as { spatial?: { coordinate_system?: string } })
            ?.spatial?.coordinate_system ?? null;

    useEffect(() => {
        if (!spatialData || !coordinateSystem) return;
        seedRenderStackFromSpatialData(config, spatialData, coordinateSystem, chart);
    }, [chart, config, coordinateSystem, spatialData]);

    const sensors = useSensors(
        useSensor(PointerSensor, {
            activationConstraint: { distance: 8 },
        }),
    );

    const insertOptions = useMemo<InsertOption[]>(() => {
        if (!stack || !spatialData) return [];
        const options: InsertOption[] = [];
        if (coordinateSystem) {
            const available = listAvailableSpatialEntries(spatialData, coordinateSystem);
            for (const entry of available) {
                if (entry.kind !== "spatial") continue;
                if (stack.entries.some((item) => item.id === entry.id)) continue;
                options.push({
                    kind: "spatial",
                    type: entry.source.elementType,
                    elementKey: entry.source.elementKey,
                    label: `${entry.source.elementType}: ${entry.source.elementKey}`,
                });
            }
        }
        for (const deckId of DECK_OVERLAY_IDS) {
            const id = deckHostLayerId(deckId);
            if (stack.entries.some((entry) => entry.id === id)) continue;
            options.push({
                kind: "host",
                deckId,
                label: DECK_OVERLAY_LABELS[deckId],
            });
        }
        return options.sort((a, b) => a.label.localeCompare(b.label));
    }, [coordinateSystem, spatialData, stack]);

    const onDragEnd = useCallback(
        (event: DragEndEvent) => {
            if (!stack) return;
            const { active, over } = event;
            if (!over || active.id === over.id) return;
            const ids = renderStackEntryIds(stack);
            const oldIndex = ids.indexOf(String(active.id));
            const newIndex = ids.indexOf(String(over.id));
            if (oldIndex === -1 || newIndex === -1) return;
            mutateStack((current) => {
                reorderRenderStackEntries(current, oldIndex, newIndex);
            });
        },
        [mutateStack, stack],
    );

    const onInsert = useCallback(
        (option: InsertOption | null) => {
            if (!option) return;
            mutateStack((current) => {
                if (option.kind === "spatial") {
                    insertSpatialRenderStackEntry(
                        current,
                        option.type,
                        option.elementKey,
                        defaultPropsForSpatialElement(option.type),
                    );
                } else {
                    insertHostRenderStackEntry(current, option.deckId);
                }
            });
        },
        [mutateStack],
    );

    if (!stack?.entries.length) {
        return (
            <div
                className="block w-full p-3"
                style={{ display: "block", alignItems: "stretch" }}
            >
                <Typography color="text.secondary">
                    Layer stack is loading…
                </Typography>
            </div>
        );
    }

    const entryIds = renderStackEntryIds(stack);

    return (
        <div
            className="block w-full space-y-2 overflow-y-auto p-3"
            style={{ display: "block", alignItems: "stretch" }}
        >
            <DndContext sensors={sensors} collisionDetection={closestCenter} onDragEnd={onDragEnd}>
                <SortableContext items={entryIds} strategy={verticalListSortingStrategy}>
                    {entryIds.map((entryId) => (
                        <SortableLayerAccordion
                            key={entryId}
                            entryId={entryId}
                            onRemove={(id) => {
                                mutateStack((current) => {
                                    removeRenderStackEntry(current, id);
                                });
                            }}
                        />
                    ))}
                </SortableContext>
            </DndContext>
            <div className="mt-4">
                <Autocomplete
                    size="small"
                    options={insertOptions}
                    groupBy={(option) =>
                        option.kind === "spatial" ? "Spatial elements" : "Deck overlays"
                    }
                    getOptionLabel={(option) => option.label}
                    onChange={(_, value) => onInsert(value)}
                    renderInput={(params) => (
                        <TextField
                            {...params}
                            label="Insert layer"
                            placeholder="Choose a layer to add"
                        />
                    )}
                />
            </div>
            <Typography variant="caption" color="text.secondary" className="mt-2 block">
                Layer order: {renderStackOrderLabel(stack)}
            </Typography>
        </div>
    );
});

export default SpatialLayerDialogComponent;
