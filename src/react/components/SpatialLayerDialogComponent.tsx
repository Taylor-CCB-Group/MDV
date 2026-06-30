import ChevronRightIcon from "@mui/icons-material/ChevronRight";
import DragIndicatorIcon from "@mui/icons-material/DragIndicator";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import HighlightOffIcon from "@mui/icons-material/HighlightOff";
import {
    Autocomplete,
    Checkbox,
    Collapse,
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
import { useEffect, useMemo, useState } from "react";
import { ErrorBoundary } from "react-error-boundary";


import {
    DECK_OVERLAY_IDS,
    DECK_OVERLAY_LABELS,
    deckHostLayerId,
    deckIdFromHostLayerId,
    type DeckOverlayId,
} from "@/react/spatialdata/host_overlay_ids";
import { renderStackEntryDisplayName } from "@/react/spatialdata/render_stack_display";
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
import { useChart, useDataStore } from "../context";
import { useConfig } from "../hooks";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "./SpatialDataMDVReact";
import DeckOverlayLayerPanel from "./spatialLayers/DeckOverlayLayerPanel";
import ImageLayerPanel from "./spatialLayers/ImageLayerPanel";
import LabelsLayerPanel from "./spatialLayers/LabelsLayerPanel";
import PointsLayerPanel from "./spatialLayers/PointsLayerPanel";
import ShapesLayerPanel from "./spatialLayers/ShapesLayerPanel";
import ErrorComponentReactWrapper from "./ErrorComponentReactWrapper";

type InsertOption =
    | { kind: "spatial"; type: RenderStackSpatialElementType; elementKey: string; label: string }
    | { kind: "host"; deckId: DeckOverlayId; label: string };

function getAvailableFields(dataStore: ReturnType<typeof useDataStore>): string[] {
    return Object.keys(dataStore.columnIndex).sort();
}

const LayerDetails = observer(function LayerDetails({ entryId }: { entryId: string }) {
    const dataStore = useDataStore();
    const chartConfig = useConfig<SpatialDataMdvReactConfig>();
    const { entry, layer, patchLayer } = useRenderStackEntry(entryId);
    const availableFields = getAvailableFields(dataStore);
    const chartColorBy =
        typeof chartConfig.color_by === "string" ? chartConfig.color_by : undefined;

    if (!entry) return null;

    if (entry.kind === "host") {
        const deckId = deckIdFromHostLayerId(entry.source.hostLayerId);
        if (!deckId) return null;
        return <DeckOverlayLayerPanel deckId={deckId} />;
    }

    if (entry.kind !== "spatial" || !layer) return null;

    switch (entry.source.elementType) {
        case "image":
            return <ImageLayerPanel entryId={entryId} />;
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
            return ( null //nothing useful here at the moment.
                // <LabelsLayerPanel
                //     config={layer as Extract<LayerConfig, { type: "labels" }>}
                //     association={NO_TABLE_ASSOCIATION}
                //     availableFields={availableFields}
                //     updateLayer={patchLayer}
                // />
            );
        default:
            return null;
    }
});

const LayerTreeRow = observer(function LayerTreeRow({
    entryId,
    selected,
    onSelect,
    onRemove,
    indent,
}: {
    entryId: string;
    selected: boolean;
    onSelect: (entryId: string) => void;
    onRemove: (entryId: string) => void;
    indent?: boolean;
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

    const layerOpacity = typeof entry.props.opacity === "number" ? entry.props.opacity : 1;
    const supportsOpacity = entry.kind === "spatial";

    const style = {
        transform: CSS.Translate.toString(transform),
        transition,
        opacity: isDragging ? 0.5 : 1,
    };

    return (
        <div
            ref={setNodeRef}
            style={style}
            className={`flex items-center gap-1 py-0.5 cursor-pointer rounded ${indent ? "pl-5 pr-1" : "px-1"} ${selected ? "bg-blue-100 dark:bg-blue-900/30" : "hover:bg-gray-100 dark:hover:bg-gray-800"}`}
            onClick={() => onSelect(entryId)}
        >
            <IconButton
                {...attributes}
                {...listeners}
                size="small"
                sx={{ cursor: "grab", flexShrink: 0 }}
                onClick={(e) => e.stopPropagation()}
            >
                <DragIndicatorIcon fontSize="small" />
            </IconButton>
            <Checkbox
                size="small"
                checked={entry.visible}
                onClick={(e) => e.stopPropagation()}
                onChange={(e) => patchEntry({ visible: e.target.checked })}
                sx={{ p: 0.25, flexShrink: 0 }}
            />
            <Typography variant="body2" className="grow truncate select-none">
                {renderStackEntryDisplayName(entry)}
            </Typography>
            {supportsOpacity && (
                <div
                    className="w-20 px-1 flex-shrink-0"
                    onClick={(e) => e.stopPropagation()}
                >
                    <Slider
                        size="small"
                        min={0}
                        max={1}
                        step={0.01}
                        value={layerOpacity}
                        onChange={(_, value) => {
                            if (typeof value === "number") patchProps({ opacity: value });
                        }}
                    />
                </div>
            )}
            {isRemovableRenderStackEntry(entry) && (
                <IconButton
                    size="small"
                    sx={{ flexShrink: 0 }}
                    onClick={(e) => {
                        e.stopPropagation();
                        onRemove(entryId);
                    }}
                    aria-label="remove layer"
                >
                    <HighlightOffIcon fontSize="small" />
                </IconButton>
            )}
        </div>
    );
});

const SpatialLayerDialogComponent = observer(() => {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const config = useConfig<SpatialDataMdvReactConfig>();
    const { spatialData } = useSpatialData();
    const stack = config.renderStack;
    const mutateStack = useRenderStackMutation();
    const [selectedEntryId, setSelectedEntryId] = useState<string | null>(null);
    // note that we'll want a proper arbitrary nested tree with groups, not the specialised state we have here
    const [spatialCollapsed, setSpatialCollapsed] = useState(false);
    // actually, just hidden in this UI for now.
    const [deckCollapsed, setDeckCollapsed] = useState(false);
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
                // multiple entries for the same element are explicitly allowed
                // (also, stack is a stable mobx reference and I think removed
                // elements weren't added back to insertOptions before)
                // now that the id is more unique rather than, determined by elementKind&Key
                // this wouldn't prevent multiple instances of the object anyway.
                // if (stack.entries.some((item) => item.id === entry.id)) continue;
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

    // Plain handlers — the React Compiler memoizes their identity.
    const makeOnDragEnd = (event: DragEndEvent) => {
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
    };

    const onInsert = (option: InsertOption | null) => {
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
    };

    const onRemove = (id: string) => {
        mutateStack((current) => removeRenderStackEntry(current, id));
        if (selectedEntryId === id) setSelectedEntryId(null);
    };

    if (!stack?.entries.length) {
        return (
            <div className="block w-full p-3">
                <Typography color="text.secondary">
                    Layer stack is loading…
                </Typography>
            </div>
        );
    }

    // crude bypass of MDV layers in dialog for now.
    const entryIds = renderStackEntryIds(stack);//.filter(id => !id.startsWith("deck:"));
    const spatialEntryIds = entryIds.filter(id => !id.startsWith("deck:"));
    const hostEntryIds = [] satisfies string[];//entryIds.filter(id => id.startsWith("deck:"));

    return (
        <div className="flex flex-col w-full h-full">
            {/* Layer tree */}
            <div className="flex-none overflow-y-auto border-b border-gray-200 dark:border-gray-700" style={{ maxHeight: "40%" }}>
                {spatialEntryIds.length > 0 && (
                    <>
                        <button
                            type="button"
                            className="flex w-full items-center gap-1 px-1 py-0.5 text-left hover:bg-gray-100 dark:hover:bg-gray-800"
                            onClick={() => setSpatialCollapsed(c => !c)}
                        >
                            {spatialCollapsed ? <ChevronRightIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}
                            <Typography variant="caption" color="text.secondary" className="font-semibold uppercase tracking-wide select-none">
                                Spatial elements
                            </Typography>
                        </button>
                        <Collapse in={!spatialCollapsed}>
                            <DndContext sensors={sensors} collisionDetection={closestCenter} onDragEnd={makeOnDragEnd}>
                                <SortableContext items={spatialEntryIds} strategy={verticalListSortingStrategy}>
                                    {spatialEntryIds.map((id) => (
                                        <LayerTreeRow
                                            key={id}
                                            entryId={id}
                                            selected={selectedEntryId === id}
                                            onSelect={setSelectedEntryId}
                                            onRemove={onRemove}
                                            indent
                                        />
                                    ))}
                                </SortableContext>
                            </DndContext>
                        </Collapse>
                    </>
                )}
                {hostEntryIds.length > 0 && (
                    <>
                        <button
                            type="button"
                            className="flex w-full items-center gap-1 px-1 py-0.5 text-left hover:bg-gray-100 dark:hover:bg-gray-800"
                            onClick={() => setDeckCollapsed(c => !c)}
                        >
                            {deckCollapsed ? <ChevronRightIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}
                            <Typography variant="caption" color="text.secondary" className="font-semibold uppercase tracking-wide select-none">
                                Deck overlays
                            </Typography>
                        </button>
                        <Collapse in={!deckCollapsed}>
                            <DndContext sensors={sensors} collisionDetection={closestCenter} onDragEnd={makeOnDragEnd}>
                                <SortableContext items={hostEntryIds} strategy={verticalListSortingStrategy}>
                                    {hostEntryIds.map((id) => (
                                        <LayerTreeRow
                                            key={id}
                                            entryId={id}
                                            selected={selectedEntryId === id}
                                            onSelect={setSelectedEntryId}
                                            onRemove={onRemove}
                                            indent
                                        />
                                    ))}
                                </SortableContext>
                            </DndContext>
                        </Collapse>
                    </>
                )}
            </div>

            {/* Properties pane */}
            <div className="flex-1 overflow-y-auto p-2 min-h-0">
                {selectedEntryId ? (
                    <ErrorBoundary
                        key={selectedEntryId}
                        FallbackComponent={({ error }) => (
                            <ErrorComponentReactWrapper
                                error={{ message: error.message, stack: error.stack }}
                                title={`Error in ${selectedEntryId}`}
                                extraMetaData={{}}
                            />
                        )}
                    >
                        <LayerDetails entryId={selectedEntryId} />
                    </ErrorBoundary>
                ) : (
                    <Typography color="text.secondary" variant="body2">
                        Select a layer to view its properties.
                    </Typography>
                )}
            </div>

            {/* Insert layer */}
            <div className="flex-none p-2 border-t border-gray-200 dark:border-gray-700">
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
        </div>
    );
});

export default SpatialLayerDialogComponent;
