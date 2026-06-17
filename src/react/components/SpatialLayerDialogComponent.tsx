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
    arrayMove,
    useSortable,
    verticalListSortingStrategy,
} from "@dnd-kit/sortable";
import { CSS } from "@dnd-kit/utilities";
import { observer } from "mobx-react-lite";
import { runInAction } from "mobx";
import { useCallback, useMemo, useState } from "react";
import { useDebouncedCallback } from "use-debounce";
import type { LayerConfig, LayerType } from "@spatialdata/vis";
import { useChart, useDataStore } from "../context";
import { useConfig } from "../hooks";
import {
    DECK_OVERLAY_IDS,
    type DeckOverlayId,
    type SpatialLayerStackConfig,
    canvasLayerOrder,
    createSpatialLayerId,
    deckStackKey,
    defaultLayerConfig,
    entryDisplayName,
    isRemovableStackEntry,
} from "@/react/spatial_layer_stack";
import { useOptionalSpatialCanvasRendererContext } from "@/react/spatial_canvas_renderer_context";
import type { SpatialDataMdvReact } from "./SpatialDataMDVReact";
import type { VivMdvReactConfig } from "./VivMDVReact";
import ImageLayerPanel from "./spatialLayers/ImageLayerPanel";
import ShapesLayerPanel from "./spatialLayers/ShapesLayerPanel";
import PointsLayerPanel from "./spatialLayers/PointsLayerPanel";
import LabelsLayerPanel from "./spatialLayers/LabelsLayerPanel";
import DeckOverlayLayerPanel from "./spatialLayers/DeckOverlayLayerPanel";

type InsertOption =
    | { kind: "spatial"; type: LayerType; elementKey: string; label: string }
    | { kind: "deck"; deckId: DeckOverlayId; label: string };

function getAvailableFields(dataStore: ReturnType<typeof useDataStore>): string[] {
    return Object.keys(dataStore.columnIndex).sort();
}

function LayerDetails({
    stackKey,
    stack,
    onUpdateSpatialLayer,
    onPatchSpatialLayer,
}: {
    stackKey: string;
    stack: SpatialLayerStackConfig;
    onUpdateSpatialLayer: (layerId: string, updates: Partial<LayerConfig>) => void;
    onPatchSpatialLayer: (layerId: string, updates: Partial<LayerConfig>) => void;
}) {
    const renderer = useOptionalSpatialCanvasRendererContext();
    const dataStore = useDataStore();
    const availableFields = useMemo(() => getAvailableFields(dataStore), [dataStore]);
    const entry = stack.entries[stackKey];
    if (!entry) return null;

    if (entry.kind === "deck") {
        return <DeckOverlayLayerPanel deckId={entry.deckId} />;
    }

    const layer = stack.spatialLayers[entry.layerId];
    if (!layer) return null;

    const association = renderer?.inferTableAssociation(layer.elementKey) ?? { status: "none" as const };
    const imageDefaults = renderer?.getImageLayerLoadedData(layer.id);
    const channelNames = useMemo(() => {
        if (!imageDefaults) return [];
        const count = Math.max(
            imageDefaults.colors?.length ?? 0,
            imageDefaults.contrastLimits?.length ?? 0,
            1,
        );
        return Array.from({ length: count }, (_, index) => `Channel ${index + 1}`);
    }, [imageDefaults]);

    switch (layer.type) {
        case "image":
            return (
                <ImageLayerPanel
                    layerId={layer.id}
                    config={layer}
                    imageSource={entry.imageSource}
                    loaderDefaults={imageDefaults}
                    channelNames={channelNames}
                    updateLayer={(updates) => onUpdateSpatialLayer(layer.id, updates)}
                    patchLayer={(updates) => onPatchSpatialLayer(layer.id, updates)}
                />
            );
        case "shapes":
            return (
                <ShapesLayerPanel
                    config={layer}
                    association={association}
                    availableFields={availableFields}
                    updateLayer={(updates) => onUpdateSpatialLayer(layer.id, updates)}
                />
            );
        case "points":
            return (
                <PointsLayerPanel
                    config={layer}
                    updateLayer={(updates) => onUpdateSpatialLayer(layer.id, updates)}
                />
            );
        case "labels":
            return (
                <LabelsLayerPanel
                    config={layer}
                    association={association}
                    availableFields={availableFields}
                    updateLayer={(updates) => onUpdateSpatialLayer(layer.id, updates)}
                />
            );
        default:
            return null;
    }
}

const SortableLayerAccordion = observer(function SortableLayerAccordion({
    stackKey,
    stack,
    onToggleVisible,
    onOpacityChange,
    onRemove,
    onUpdateSpatialLayer,
    onPatchSpatialLayer,
}: {
    stackKey: string;
    stack: SpatialLayerStackConfig;
    onToggleVisible: (key: string, visible: boolean) => void;
    onOpacityChange: (key: string, opacity: number) => void;
    onRemove: (key: string) => void;
    onUpdateSpatialLayer: (layerId: string, updates: Partial<LayerConfig>) => void;
    onPatchSpatialLayer: (layerId: string, updates: Partial<LayerConfig>) => void;
}) {
    const {
        attributes,
        listeners,
        setNodeRef,
        transform,
        transition,
        isDragging,
    } = useSortable({ id: stackKey });
    const [isHovered, setIsHovered] = useState(false);
    const entry = stack.entries[stackKey];
    if (!entry) return null;

    const opacity = entry.opacity ?? 1;
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
                {isRemovableStackEntry(stackKey, stack) && (
                    <IconButton
                        onClick={(event) => {
                            event.stopPropagation();
                            onRemove(stackKey);
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
                        onChange={(event) => onToggleVisible(stackKey, event.target.checked)}
                    />
                    <Typography variant="subtitle1" className="grow">
                        {entryDisplayName(stackKey, stack)}
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
                                        onOpacityChange(stackKey, value);
                                    }
                                }}
                            />
                        </div>
                    )}
                </div>
            </AccordionSummary>
            <AccordionDetails>
                <LayerDetails
                    stackKey={stackKey}
                    stack={stack}
                    onUpdateSpatialLayer={onUpdateSpatialLayer}
                    onPatchSpatialLayer={onPatchSpatialLayer}
                />
            </AccordionDetails>
        </Accordion>
    );
});

const SpatialLayerDialogComponent = observer(() => {
    const chart = useChart<VivMdvReactConfig, SpatialDataMdvReact>();
    const config = useConfig<VivMdvReactConfig & { spatialLayerStack?: SpatialLayerStackConfig }>();
    const renderer = useOptionalSpatialCanvasRendererContext();
    const stack = config.spatialLayerStack;
    if (!stack?.stackOrder.length) {
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

    const sensors = useSensors(
        useSensor(PointerSensor, {
            activationConstraint: { distance: 8 },
        }),
    );

    const updateStack = useCallback(
        (mutator: (current: SpatialLayerStackConfig) => void, syncCanvas = false) => {
            runInAction(() => {
                const current = config.spatialLayerStack;
                if (!current) return;
                mutator(current);
                if (syncCanvas) {
                    chart.updateCanvasLayerState(current);
                }
            });
        },
        [chart, config],
    );

    const persistStack = useDebouncedCallback(
        (mutator: (current: SpatialLayerStackConfig) => void) => {
            updateStack(mutator);
        },
        150,
    );

    const patchSpatialLayer = useCallback(
        (layerId: string, updates: Partial<LayerConfig>) => {
            chart.patchCanvasLayer(layerId, updates);
        },
        [chart],
    );

    const persistSpatialLayer = useCallback(
        (layerId: string, updates: Partial<LayerConfig>) => {
            persistStack((current) => {
                const layer = current.spatialLayers[layerId];
                if (!layer) return;
                if ("channels" in updates && updates.channels && "channels" in layer) {
                    layer.channels = {
                        ...layer.channels,
                        ...updates.channels,
                    };
                }
                const rest = { ...updates };
                if ("channels" in rest) {
                    delete rest.channels;
                }
                Object.assign(layer, rest);
            });
        },
        [persistStack],
    );

    const insertOptions = useMemo<InsertOption[]>(() => {
        const options: InsertOption[] = [];
        const available = renderer?.availableElements;
        if (available) {
            const buckets: Array<{ type: LayerType; key: keyof typeof available }> = [
                { type: "image", key: "images" },
                { type: "shapes", key: "shapes" },
                { type: "points", key: "points" },
                { type: "labels", key: "labels" },
            ];
            for (const { type, key } of buckets) {
                for (const element of available[key]) {
                    const layerId = createSpatialLayerId(type, element.key);
                    if (stack.spatialLayers[layerId]) continue;
                    options.push({
                        kind: "spatial",
                        type,
                        elementKey: element.key,
                        label: `${type}: ${element.key}`,
                    });
                }
            }
        }
        for (const deckId of DECK_OVERLAY_IDS) {
            const key = deckStackKey(deckId);
            if (stack.entries[key]) continue;
            options.push({
                kind: "deck",
                deckId,
                label: `deck: ${deckId}`,
            });
        }
        return options.sort((a, b) => a.label.localeCompare(b.label));
    }, [renderer?.availableElements, stack.entries, stack.spatialLayers]);

    const onDragEnd = useCallback(
        (event: DragEndEvent) => {
            const { active, over } = event;
            if (!over || active.id === over.id) return;
            const oldIndex = stack.stackOrder.indexOf(String(active.id));
            const newIndex = stack.stackOrder.indexOf(String(over.id));
            if (oldIndex === -1 || newIndex === -1) return;
            updateStack((current) => {
                current.stackOrder.splice(
                    0,
                    current.stackOrder.length,
                    ...arrayMove([...current.stackOrder], oldIndex, newIndex),
                );
            }, true);
        },
        [stack.stackOrder, updateStack],
    );

    const onInsert = useCallback(
        (option: InsertOption | null) => {
            if (!option) return;
            updateStack((current) => {
                if (option.kind === "spatial") {
                    const layerId = createSpatialLayerId(option.type, option.elementKey);
                    current.spatialLayers[layerId] =
                        current.spatialLayers[layerId] ??
                        defaultLayerConfig(option.type, option.elementKey);
                    current.entries[layerId] = {
                        kind: "spatial",
                        layerId,
                        visible: true,
                        opacity: 1,
                        imageSource: option.type === "image" ? "spatial" : undefined,
                    };
                    if (!current.stackOrder.includes(layerId)) {
                        current.stackOrder.push(layerId);
                    }
                } else {
                    const key = deckStackKey(option.deckId);
                    current.entries[key] = {
                        kind: "deck",
                        deckId: option.deckId,
                        visible: true,
                        opacity: 1,
                    };
                    if (!current.stackOrder.includes(key)) {
                        current.stackOrder.push(key);
                    }
                }
            }, true);
        },
        [updateStack],
    );

    return (
        <div
            className="block w-full space-y-2 overflow-y-auto p-3"
            style={{ display: "block", alignItems: "stretch" }}
        >
            <DndContext sensors={sensors} collisionDetection={closestCenter} onDragEnd={onDragEnd}>
                <SortableContext
                    items={stack.stackOrder}
                    strategy={verticalListSortingStrategy}
                >
                    {stack.stackOrder.map((stackKey: string) => (
                        <SortableLayerAccordion
                            key={stackKey}
                            stackKey={stackKey}
                            stack={stack}
                            onToggleVisible={(key, visible) => {
                                const entry = stack.entries[key];
                                if (entry?.kind === "spatial") {
                                    patchSpatialLayer(entry.layerId, { visible });
                                }
                                updateStack((current) => {
                                    const currentEntry = current.entries[key];
                                    if (!currentEntry) return;
                                    currentEntry.visible = visible;
                                }, entry?.kind === "deck");
                            }}
                            onOpacityChange={(key, opacity) => {
                                const entry = stack.entries[key];
                                if (entry?.kind === "spatial") {
                                    patchSpatialLayer(entry.layerId, { opacity });
                                    persistStack((current) => {
                                        const currentEntry = current.entries[key];
                                        if (!currentEntry) return;
                                        currentEntry.opacity = opacity;
                                        const layer = current.spatialLayers[entry.layerId];
                                        if (layer) {
                                            layer.opacity = opacity;
                                        }
                                    });
                                }
                            }}
                            onRemove={(key) => {
                                updateStack((current) => {
                                    const entry = current.entries[key];
                                    delete current.entries[key];
                                    if (entry?.kind === "spatial") {
                                        delete current.spatialLayers[entry.layerId];
                                    }
                                    const nextOrder = current.stackOrder.filter((id) => id !== key);
                                    current.stackOrder.splice(
                                        0,
                                        current.stackOrder.length,
                                        ...nextOrder,
                                    );
                                }, true);
                            }}
                            onUpdateSpatialLayer={(layerId, updates) => {
                                persistSpatialLayer(layerId, updates);
                            }}
                            onPatchSpatialLayer={patchSpatialLayer}
                        />
                    ))}
                </SortableContext>
            </DndContext>
            <div className="mt-4">
                <Autocomplete
                    size="small"
                    options={insertOptions}
                    groupBy={(option) => (option.kind === "spatial" ? "Spatial elements" : "Deck overlays")}
                    getOptionLabel={(option) => option.label}
                    onChange={(_, value) => onInsert(value)}
                    renderInput={(params) => (
                        <TextField {...params} label="Insert layer" placeholder="Choose a layer to add" />
                    )}
                />
            </div>
            <Typography variant="caption" color="text.secondary" className="mt-2 block">
                Layer order: {canvasLayerOrder(stack).join(" → ")}
            </Typography>
        </div>
    );
});

export default SpatialLayerDialogComponent;
