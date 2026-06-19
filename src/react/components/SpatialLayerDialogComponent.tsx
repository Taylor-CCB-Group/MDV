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
import type { RenderStack, RenderStackEntry, RenderStackSpatialElementType } from "@spatialdata/layers";
import type { LayerConfig } from "@spatialdata/vis";
import { useSpatialData } from "@spatialdata/react";
import { observer } from "mobx-react-lite";
import { runInAction } from "mobx";
import { useCallback, useMemo, useState } from "react";

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
    patchRenderStackEntry,
    removeRenderStackEntry,
    reorderRenderStackEntries,
    renderStackEntryIds,
} from "@/react/spatialdata/render_stack_mutations";
import {
    defaultPropsForSpatialElement,
    listAvailableSpatialEntries,
} from "@/react/spatialdata/render_stack_seed";
import { NO_TABLE_ASSOCIATION } from "@/react/spatialdata/table_association";
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

function spatialPropsAsLayerConfig(
    entry: Extract<RenderStackEntry, { kind: "spatial" }>,
): LayerConfig {
    return {
        id: entry.id,
        type: entry.source.elementType,
        elementKey: entry.source.elementKey,
        visible: entry.visible,
        opacity: typeof entry.props.opacity === "number" ? entry.props.opacity : 1,
        ...entry.props,
    } as LayerConfig;
}

function LayerDetails({
    entry,
    onPatchProps,
}: {
    entry: RenderStackEntry;
    onPatchProps: (entryId: string, props: Record<string, unknown>, merge?: boolean) => void;
}) {
    const dataStore = useDataStore();
    const availableFields = useMemo(() => getAvailableFields(dataStore), [dataStore]);

    if (entry.kind === "host") {
        const deckId = deckIdFromHostLayerId(entry.source.hostLayerId);
        if (!deckId) return null;
        return <DeckOverlayLayerPanel deckId={deckId} />;
    }

    if (entry.kind !== "spatial") return null;

    const layer = spatialPropsAsLayerConfig(entry);
    const patchLayer = (updates: Partial<LayerConfig>) => {
        onPatchProps(entry.id, updates as Record<string, unknown>);
    };

    switch (entry.source.elementType) {
        case "image":
            return (
                <ImageLayerPanel
                    layerId={entry.id}
                    config={layer as Extract<LayerConfig, { type: "image" }>}
                    updateLayer={patchLayer}
                    patchLayer={patchLayer}
                />
            );
        case "shapes":
            return (
                <ShapesLayerPanel
                    config={layer as Extract<LayerConfig, { type: "shapes" }>}
                    association={NO_TABLE_ASSOCIATION}
                    availableFields={availableFields}
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
}

const SortableLayerAccordion = observer(function SortableLayerAccordion({
    entry,
    onToggleVisible,
    onOpacityChange,
    onRemove,
    onPatchProps,
}: {
    entry: RenderStackEntry;
    onToggleVisible: (entryId: string, visible: boolean) => void;
    onOpacityChange: (entryId: string, opacity: number) => void;
    onRemove: (entryId: string) => void;
    onPatchProps: (entryId: string, props: Record<string, unknown>) => void;
}) {
    const {
        attributes,
        listeners,
        setNodeRef,
        transform,
        transition,
        isDragging,
    } = useSortable({ id: entry.id });
    const [isHovered, setIsHovered] = useState(false);

    const opacity =
        typeof entry.props.opacity === "number" ? entry.props.opacity : 1;
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
            onMouseEnter={() => setIsHovered(true)}
            onMouseLeave={() => setIsHovered(false)}
        >
            <AccordionSummary expandIcon={<ArrowDropDownIcon />}>
                {isRemovableRenderStackEntry(entry) && (
                    <IconButton
                        onClick={(event) => {
                            event.stopPropagation();
                            onRemove(entry.id);
                        }}
                        aria-label="remove layer"
                        size="small"
                        sx={{
                            position: "absolute",
                            right: "-18px",
                            top: "-18px",
                            opacity: isHovered ? 1 : 0,
                            transition: "opacity 0.3s",
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
                        onChange={(event) => onToggleVisible(entry.id, event.target.checked)}
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
                                        onOpacityChange(entry.id, value);
                                    }
                                }}
                            />
                        </div>
                    )}
                </div>
            </AccordionSummary>
            <AccordionDetails>
                <LayerDetails entry={entry} onPatchProps={onPatchProps} />
            </AccordionDetails>
        </Accordion>
    );
});

const SpatialLayerDialogComponent = observer(() => {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const config = useConfig<SpatialDataMdvReactConfig>();
    const { spatialData } = useSpatialData();
    const stack = config.renderStack;

    const sensors = useSensors(
        useSensor(PointerSensor, {
            activationConstraint: { distance: 8 },
        }),
    );

    const updateStack = useCallback(
        (mutator: (current: RenderStack) => void) => {
            runInAction(() => {
                if (!config.renderStack) return;
                mutator(config.renderStack);
            });
        },
        [config],
    );

    const onPatchProps = useCallback(
        (entryId: string, props: Record<string, unknown>) => {
            updateStack((current) => {
                patchRenderStackEntry(current, entryId, { props });
            });
        },
        [updateStack],
    );

    const insertOptions = useMemo<InsertOption[]>(() => {
        if (!stack || !spatialData) return [];
        const options: InsertOption[] = [];
        const coordinateSystem =
            (chart.dataStore.regions?.all_regions[config.region] as { spatial?: { coordinate_system?: string } })
                ?.spatial?.coordinate_system ?? null;
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
    }, [chart.dataStore.regions?.all_regions, config.region, spatialData, stack]);

    const onDragEnd = useCallback(
        (event: DragEndEvent) => {
            if (!stack) return;
            const { active, over } = event;
            if (!over || active.id === over.id) return;
            const ids = renderStackEntryIds(stack);
            const oldIndex = ids.indexOf(String(active.id));
            const newIndex = ids.indexOf(String(over.id));
            if (oldIndex === -1 || newIndex === -1) return;
            updateStack((current) => {
                reorderRenderStackEntries(current, oldIndex, newIndex);
            });
        },
        [stack, updateStack],
    );

    const onInsert = useCallback(
        (option: InsertOption | null) => {
            if (!option) return;
            updateStack((current) => {
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
        [updateStack],
    );

    if (!stack?.entries.length) {
        return (
            <div
                className="block w-full p-3"
                style={{ display: "block", alignItems: "stretch" }}
            >
                <Typography color="text.secondary">
                    Layer stack is loading. Open this dialog again after spatial data has loaded.
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
                    {stack.entries.map((entry) => (
                        <SortableLayerAccordion
                            key={entry.id}
                            entry={entry}
                            onToggleVisible={(entryId, visible) => {
                                updateStack((current) => {
                                    patchRenderStackEntry(current, entryId, { visible });
                                });
                            }}
                            onOpacityChange={(entryId, opacity) => {
                                updateStack((current) => {
                                    patchRenderStackEntry(current, entryId, {
                                        props: { opacity },
                                    });
                                });
                            }}
                            onRemove={(entryId) => {
                                updateStack((current) => {
                                    removeRenderStackEntry(current, entryId);
                                });
                            }}
                            onPatchProps={onPatchProps}
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
