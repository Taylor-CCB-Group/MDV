import {
    closestCenter,
    DndContext,
    type DragEndEvent,
    type DragMoveEvent,
    type DragOverEvent,
    type DragStartEvent,
    PointerSensor,
    useSensor,
    useSensors,
} from "@dnd-kit/core";
import { SortableContext, useSortable, verticalListSortingStrategy } from "@dnd-kit/sortable";
import { CSS } from "@dnd-kit/utilities";
import AddIcon from "@mui/icons-material/Add";
import ChevronRightIcon from "@mui/icons-material/ChevronRight";
import DragIndicatorIcon from "@mui/icons-material/DragIndicator";
import ErrorOutlineIcon from "@mui/icons-material/ErrorOutline";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import HighlightOffIcon from "@mui/icons-material/HighlightOff";
import InfoOutlinedIcon from "@mui/icons-material/InfoOutlined";
import LayersClearIcon from "@mui/icons-material/LayersClear";
import WarningAmberIcon from "@mui/icons-material/WarningAmber";
import {
    Alert,
    AlertTitle,
    Checkbox,
    Collapse,
    Divider,
    IconButton,
    ListSubheader,
    Menu,
    MenuItem,
    Slider,
    TextField,
    Typography,
} from "@mui/material";
import { findRenderStackEntry, type RenderStackSpatialElementType, walkRenderStackEntries } from "@spatialdata/layers";
import type { LayerConfig } from "@spatialdata/vis";
import { observer } from "mobx-react-lite";
import { type KeyboardEvent, useCallback, useEffect, useMemo, useState } from "react";
import { ErrorBoundary } from "react-error-boundary";
import { compositingBadge } from "@/react/spatialdata/blend_modes";
import {
    DECK_OVERLAY_VISIBILITY_FOLLOWERS,
    type DeckOverlayId,
    deckHostLayerId,
    deckIdFromHostLayerId,
    isDeckHostLayerId,
    LISTED_DECK_OVERLAYS,
} from "@/react/spatialdata/host_overlay_ids";
import { flattenLayerTree, projectLayerTreeDrop } from "@/react/spatialdata/layer_tree";
import {
    findRenderStackSiblings,
    insertRenderStackGroup,
    insertSpatialRenderStackEntry,
    isRemovableRenderStackEntry,
    moveRenderStackEntry,
    patchRenderStackEntry,
    removeRenderStackEntry,
    seedRenderStackFromSpatialData,
    ungroupRenderStackEntry,
    useRenderStackEntry,
    useRenderStackMutation,
} from "@/react/spatialdata/render_stack_control";
import { defaultPropsForSpatialElement, listAvailableSpatialEntries } from "@/react/spatialdata/render_stack_defaults";
import {
    customDisplayName,
    DISPLAY_NAME_PROP,
    numberedDefaultName,
    renderStackEntryDisplayName,
} from "@/react/spatialdata/render_stack_display";
import { useChartSpatialData } from "@/react/spatialdata/shared_spatial_data";
import { useElementTableAssociation } from "@/react/spatialdata/table_association";
import { useChart, useDataStore } from "../context";
import { useConfig, useDataSources } from "../hooks";
import EntryStatusDetails from "./EntryStatusDetails";
import ErrorComponentReactWrapper from "./ErrorComponentReactWrapper";
import type { SpatialDataMdvReact, SpatialDataMdvReactConfig } from "./SpatialDataMDVReact";
import DeckOverlayLayerPanel from "./spatialLayers/DeckOverlayLayerPanel";
import GroupLayerPanel from "./spatialLayers/GroupLayerPanel";
import ImageLayerPanel from "./spatialLayers/ImageLayerPanel";
import LabelsLayerPanel from "./spatialLayers/LabelsLayerPanel";
import LayerBlendModePanel from "./spatialLayers/LayerBlendModePanel";
import PointsLayerPanel from "./spatialLayers/PointsLayerPanel";
import ShapesLayerPanel from "./spatialLayers/ShapesLayerPanel";

/** SpatialData elements, or an empty group: MDV overlays are always in the stack, shown or hidden. */
type InsertOption =
    | { kind: "spatial"; type: RenderStackSpatialElementType; elementKey: string; label: string }
    | { kind: "group"; label: string };

const GROUP_INSERT_OPTION: InsertOption = { kind: "group", label: "Group" };

/** How far each level of the layer tree is indented, and how far a drag must go sideways to change level. */
const TREE_INDENT_PX = 20;

function isLayerConfigOfType<TType extends LayerConfig["type"]>(
    layer: LayerConfig | undefined,
    type: TType,
): layer is Extract<LayerConfig, { type: TType }> {
    return layer?.type === type;
}

/**
 * What sd.js reports about the layer: why it could not load (with a retry where that could
 * help), and what is degraded but still drawing.
 */
const LayerStatusDetails = observer(function LayerStatusDetails({ entryId }: { entryId: string }) {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const status = chart.entryStatuses.get(entryId);
    if (!status) return null;
    const { retryEntry } = chart;
    return (
        <Alert
            severity={status.errors.length > 0 ? "error" : "info"}
            sx={{ mb: 1, "& .MuiAlert-message": { minWidth: 0, width: "100%" } }}
        >
            <AlertTitle>{status.errors.length > 0 ? "This layer could not load" : "Notes on this layer"}</AlertTitle>
            <Typography variant="body2" component="div">
                <EntryStatusDetails status={status} onRetry={retryEntry && (() => void retryEntry(entryId))} />
            </Typography>
        </Alert>
    );
});

/** Why a layer was turned off; shown above its properties until it is turned back on. */
const LayerFailureDetails = observer(function LayerFailureDetails({ entryId }: { entryId: string }) {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const failure = chart.layerFailures.get(entryId);
    if (!failure) return null;
    return (
        <Alert severity="warning" sx={{ mb: 1, "& .MuiAlert-message": { minWidth: 0 } }}>
            <AlertTitle>This layer failed to draw and was turned off</AlertTitle>
            <Typography variant="body2" sx={{ wordBreak: "break-word" }}>
                {failure.reason}
            </Typography>
            {failure.detail && (
                <Typography
                    variant="caption"
                    component="pre"
                    sx={{ mt: 1, whiteSpace: "pre-wrap", wordBreak: "break-word", opacity: 0.8 }}
                >
                    {failure.detail}
                </Typography>
            )}
            <Typography variant="caption" component="p" sx={{ mt: 1 }}>
                {failure.source === "webgpu" ? "WebGPU" : "deck.gl"} error
                {failure.count > 1 ? `, reported ${failure.count} times` : ""}. Turn the layer back on to try again.
            </Typography>
        </Alert>
    );
});

const LayerDetails = observer(function LayerDetails({ entryId }: { entryId: string }) {
    const dataStore = useDataStore();
    const dataSources = useDataSources();
    const { spatialData } = useChartSpatialData();
    const { entry, layer, patchLayer } = useRenderStackEntry(entryId);
    const spatialElementType = entry?.kind === "spatial" ? entry.source.elementType : undefined;
    const spatialElementKey = entry?.kind === "spatial" ? entry.source.elementKey : undefined;
    const tableAssociation = useElementTableAssociation(
        spatialData ?? undefined,
        spatialElementType,
        spatialElementKey,
        dataSources,
    );
    const layerDataStore =
        tableAssociation.status === "resolved"
            ? (dataSources.find((dataSource) => dataSource.name === tableAssociation.dataSourceName)?.dataStore ??
              dataStore)
            : dataStore;

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
            if (!isLayerConfigOfType(layer, "shapes")) return null;
            return (
                <ShapesLayerPanel
                    config={layer}
                    association={tableAssociation}
                    dataStore={layerDataStore}
                    updateLayer={patchLayer}
                />
            );
        case "points":
            if (!isLayerConfigOfType(layer, "points")) return null;
            return <PointsLayerPanel config={layer} updateLayer={patchLayer} />;
        case "labels":
            if (!isLayerConfigOfType(layer, "labels")) return null;
            return (
                <LabelsLayerPanel
                    config={layer}
                    association={tableAssociation}
                    dataStore={layerDataStore}
                    updateLayer={patchLayer}
                />
            );
        default:
            return null;
    }
});

/**
 * Edits an entry's name. The draft is local and committed on Enter or blur, so typing
 * does not patch the stack once per keystroke; an empty name goes back to the default.
 */
const LayerNameField = observer(function LayerNameField({
    entryId,
    autoFocus,
    onDone,
    variant = "outlined",
}: {
    entryId: string;
    autoFocus?: boolean;
    onDone?: () => void;
    variant?: "outlined" | "standard";
}) {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const { entry, patchProps } = useRenderStackEntry(entryId);
    const saved = entry ? (customDisplayName(entry) ?? "") : "";
    const [draft, setDraft] = useState(saved);
    // follow renames made elsewhere (the other name field) while this one is not being edited
    const [editing, setEditing] = useState(false);
    useEffect(() => {
        if (!editing) setDraft(saved);
    }, [saved, editing]);
    if (!entry) return null;
    const placeholder = numberedDefaultName(entry, config.renderStack);
    const commit = () => {
        setEditing(false);
        const next = draft.trim();
        if (next !== saved) patchProps({ [DISPLAY_NAME_PROP]: next === "" ? undefined : next });
        onDone?.();
    };
    const onKeyDown = (e: KeyboardEvent) => {
        e.stopPropagation();
        // Enter that confirms an IME composition (e.g. Japanese input) is not the end of the name
        if (e.key === "Enter" && !e.nativeEvent.isComposing) commit();
        if (e.key === "Escape") {
            setEditing(false);
            setDraft(saved);
            onDone?.();
        }
    };
    return (
        <TextField
            size="small"
            fullWidth
            variant={variant}
            label={variant === "outlined" ? "Name" : undefined}
            autoFocus={autoFocus}
            value={draft}
            placeholder={placeholder}
            onFocus={() => setEditing(true)}
            onChange={(e) => {
                setEditing(true);
                setDraft(e.target.value);
            }}
            onBlur={commit}
            onKeyDown={onKeyDown}
            onClick={(e) => e.stopPropagation()}
            slotProps={{ htmlInput: { "aria-label": "layer name" } }}
        />
    );
});

/** What the selected layer is, above its properties: its name, and what it draws. */
const LayerHeader = observer(function LayerHeader({ entryId }: { entryId: string }) {
    const { entry } = useRenderStackEntry(entryId);
    if (!entry) return null;
    const source =
        entry.kind === "spatial"
            ? `SpatialData ${entry.source.elementType} element '${entry.source.elementKey}'`
            : entry.kind === "host"
              ? "MDV overlay"
              : "Layer group";
    return (
        <div className="grid gap-1 mb-3">
            <LayerNameField entryId={entryId} />
            {source && (
                <Typography variant="caption" color="text.secondary">
                    {source}
                </Typography>
            )}
        </div>
    );
});

/** Ids whose visibility a row's checkbox sets: its own, plus any overlay that follows it. */
function visibilityTargets(entryId: string): string[] {
    const deckId = deckIdFromHostLayerId(entryId);
    const followers = deckId ? (DECK_OVERLAY_VISIBILITY_FOLLOWERS[deckId] ?? []) : [];
    return [entryId, ...followers.map(deckHostLayerId)];
}

/** A visibility checkbox's change: the entries (and their followers) shown or hidden together. */
function useSetEntriesVisible() {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const mutateStack = useRenderStackMutation();
    return useCallback(
        (entryIds: string[], visible: boolean) => {
            const targets = entryIds.flatMap(visibilityTargets);
            // turning a failed layer back on is a retry: forget the failure so it can be reported again
            if (visible) for (const id of targets) chart.clearLayerFailure(id);
            mutateStack((stack) => {
                for (const id of targets) patchRenderStackEntry(stack, id, { visible });
            });
        },
        [chart, mutateStack],
    );
}

type DragHandleProps = Pick<ReturnType<typeof useSortable>, "attributes" | "listeners">;

const LayerTreeRow = observer(function LayerTreeRow({
    entryId,
    selected,
    onSelect,
    onRemove,
    dragHandle,
    setNodeRef,
    style,
    nested = false,
    depth = 0,
    expander,
    onUngroup,
    onAddInto,
}: {
    entryId: string;
    selected: boolean;
    onSelect: (entryId: string) => void;
    onRemove: (entryId: string) => void;
    /** Absent for rows that cannot be reordered. */
    dragHandle?: DragHandleProps;
    setNodeRef?: (node: HTMLElement | null) => void;
    style?: React.CSSProperties;
    /** Inside a dialog-only overlay group row. */
    nested?: boolean;
    /** How deep in the stack's groups the row sits. */
    depth?: number;
    /** For a layer group's row: whether its children are shown. */
    expander?: { expanded: boolean; onToggle: () => void };
    onUngroup?: (entryId: string) => void;
    /** For a layer group's row: open the add-layer menu, adding into this group. */
    onAddInto?: (groupId: string, anchor: HTMLElement) => void;
}) {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const config = useConfig<SpatialDataMdvReactConfig>();
    const setEntriesVisible = useSetEntriesVisible();
    const { entry, patchProps } = useRenderStackEntry(entryId);
    const [renaming, setRenaming] = useState(false);
    const failure = chart.layerFailures.get(entryId);
    const status = chart.entryStatuses.get(entryId);

    if (!entry) return null;

    const layerOpacity = typeof entry.props.opacity === "number" ? entry.props.opacity : 1;
    const supportsOpacity = entry.kind === "spatial" || entry.kind === "group";
    const name = renderStackEntryDisplayName(entry, config.renderStack);
    const badge = config.renderStack && compositingBadge(config.renderStack, entry);

    return (
        <div
            ref={setNodeRef}
            style={{ paddingLeft: (nested ? 40 : 20) + depth * TREE_INDENT_PX, ...style }}
            className={`flex items-center gap-1 py-0.5 cursor-pointer rounded pr-1 ${selected ? "bg-blue-100 dark:bg-blue-900/30" : failure ? "bg-amber-100 hover:bg-amber-200 dark:bg-amber-900/30 dark:hover:bg-amber-900/50" : "hover:bg-gray-100 dark:hover:bg-gray-800"}`}
            onClick={() => onSelect(entryId)}
            title={failure ? `Failed and turned off: ${failure.reason}` : undefined}
        >
            {dragHandle ? (
                <IconButton
                    {...dragHandle.attributes}
                    {...dragHandle.listeners}
                    size="small"
                    sx={{ cursor: "grab", flexShrink: 0 }}
                    onClick={(e) => e.stopPropagation()}
                >
                    <DragIndicatorIcon fontSize="small" />
                </IconButton>
            ) : (
                <span className="w-[30px] flex-shrink-0" />
            )}
            {expander && (
                <IconButton
                    size="small"
                    sx={{ flexShrink: 0, p: 0 }}
                    onClick={(e) => {
                        e.stopPropagation();
                        expander.onToggle();
                    }}
                    aria-label={expander.expanded ? `collapse ${name}` : `expand ${name}`}
                    aria-expanded={expander.expanded}
                >
                    {expander.expanded ? <ExpandMoreIcon fontSize="small" /> : <ChevronRightIcon fontSize="small" />}
                </IconButton>
            )}
            <Checkbox
                size="small"
                checked={entry.visible}
                onClick={(e) => e.stopPropagation()}
                onChange={(e) => setEntriesVisible([entryId], e.target.checked)}
                sx={{ p: 0.25, flexShrink: 0 }}
                inputProps={{ "aria-label": `show ${name}` }}
            />
            {failure && <WarningAmberIcon fontSize="small" color="warning" sx={{ flexShrink: 0 }} />}
            {!failure && status && (
                <span
                    className="flex flex-shrink-0"
                    title={
                        status.errors.length > 0
                            ? `Could not load: ${status.errors.map(({ error }) => error.message).join(" ")}`
                            : status.notices.map((notice) => notice.message).join(" ")
                    }
                >
                    {status.errors.length > 0 ? (
                        <ErrorOutlineIcon fontSize="small" color="error" />
                    ) : (
                        <InfoOutlinedIcon fontSize="small" color="info" />
                    )}
                </span>
            )}
            {renaming ? (
                <div className="grow min-w-0">
                    <LayerNameField entryId={entryId} autoFocus variant="standard" onDone={() => setRenaming(false)} />
                </div>
            ) : (
                <Typography
                    variant="body2"
                    className="grow truncate select-none rounded focus-visible:outline focus-visible:outline-2"
                    title="Double-click (or F2) to rename"
                    // the keyboard path to selecting a row; the row's own click is for the mouse
                    role="button"
                    tabIndex={0}
                    aria-pressed={selected}
                    onKeyDown={(e) => {
                        if (e.key === "Enter" || e.key === " ") {
                            e.preventDefault();
                            onSelect(entryId);
                        } else if (e.key === "F2") {
                            e.preventDefault();
                            onSelect(entryId);
                            setRenaming(true);
                        }
                    }}
                    onDoubleClick={(e) => {
                        e.stopPropagation();
                        onSelect(entryId);
                        setRenaming(true);
                    }}
                >
                    {name}
                </Typography>
            )}
            {badge && (
                <Typography
                    variant="caption"
                    className={`flex-shrink-0 rounded border px-1 leading-4 select-none ${badge.inactive ? "border-dashed opacity-60" : "border-current"}`}
                    color={badge.inactive ? "text.secondary" : "primary"}
                    title={badge.title}
                >
                    {badge.label}
                </Typography>
            )}
            {entry.kind === "group" && onAddInto && (
                <IconButton
                    size="small"
                    sx={{ flexShrink: 0 }}
                    onClick={(e) => {
                        e.stopPropagation();
                        onAddInto(entryId, e.currentTarget);
                    }}
                    aria-label={`add a layer to ${name}`}
                    title="Add a layer to this group"
                >
                    <AddIcon fontSize="small" />
                </IconButton>
            )}
            {supportsOpacity && (
                <div className="w-20 px-1 flex-shrink-0" onClick={(e) => e.stopPropagation()}>
                    <Slider
                        size="small"
                        min={0}
                        max={1}
                        step={0.01}
                        value={layerOpacity}
                        aria-label={`${name} opacity`}
                        onChange={(_, value) => {
                            if (typeof value === "number") patchProps({ opacity: value });
                        }}
                    />
                </div>
            )}
            {entry.kind === "group" && onUngroup && (
                <IconButton
                    size="small"
                    sx={{ flexShrink: 0 }}
                    onClick={(e) => {
                        e.stopPropagation();
                        onUngroup(entryId);
                    }}
                    aria-label="ungroup"
                    title="Ungroup: keep the layers, remove the group"
                >
                    <LayersClearIcon fontSize="small" />
                </IconButton>
            )}
            {entry.kind === "spatial" && isRemovableRenderStackEntry(entry) && (
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

/**
 * A dialog-only group of MDV overlays (the stack keeps them flat): one checkbox for all
 * of them, and the overlays themselves beneath it.
 */
const OverlayGroupRow = observer(function OverlayGroupRow({
    label,
    deckIds,
    selectedEntryId,
    onSelect,
    onRemove,
}: {
    label: string;
    deckIds: DeckOverlayId[];
    selectedEntryId: string | null;
    onSelect: (entryId: string) => void;
    onRemove: (entryId: string) => void;
}) {
    const config = useConfig<SpatialDataMdvReactConfig>();
    const setEntriesVisible = useSetEntriesVisible();
    const [collapsed, setCollapsed] = useState(true);
    const stack = config.renderStack;
    const entryIds = deckIds.map(deckHostLayerId).filter((id) => stack && findRenderStackEntry(stack, id));
    if (entryIds.length === 0) return null;
    const shown = entryIds.filter((id) => stack && findRenderStackEntry(stack, id)?.visible);

    return (
        <>
            <div className="flex items-center gap-1 py-0.5 rounded pl-5 pr-1">
                <IconButton
                    size="small"
                    sx={{ flexShrink: 0 }}
                    onClick={() => setCollapsed((c) => !c)}
                    aria-label={collapsed ? `expand ${label}` : `collapse ${label}`}
                    aria-expanded={!collapsed}
                >
                    {collapsed ? <ChevronRightIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}
                </IconButton>
                <Checkbox
                    size="small"
                    checked={shown.length === entryIds.length}
                    indeterminate={shown.length > 0 && shown.length < entryIds.length}
                    onChange={(e) => setEntriesVisible(entryIds, e.target.checked)}
                    sx={{ p: 0.25, flexShrink: 0 }}
                    inputProps={{ "aria-label": `show ${label}` }}
                />
                <Typography variant="body2" className="grow truncate select-none">
                    {label}
                </Typography>
            </div>
            <Collapse in={!collapsed}>
                {entryIds.map((id) => (
                    <LayerTreeRow
                        key={id}
                        entryId={id}
                        selected={selectedEntryId === id}
                        onSelect={onSelect}
                        onRemove={onRemove}
                        nested
                    />
                ))}
            </Collapse>
        </>
    );
});

function SortableLayerTreeRow(props: {
    entryId: string;
    selected: boolean;
    onSelect: (entryId: string) => void;
    onRemove: (entryId: string) => void;
    depth: number;
    expander?: { expanded: boolean; onToggle: () => void };
    onUngroup?: (entryId: string) => void;
    onAddInto?: (groupId: string, anchor: HTMLElement) => void;
}) {
    const { attributes, listeners, setNodeRef, transform, transition, isDragging } = useSortable({
        id: props.entryId,
    });
    return (
        <LayerTreeRow
            {...props}
            dragHandle={{ attributes, listeners }}
            setNodeRef={setNodeRef}
            style={{ transform: CSS.Translate.toString(transform), transition, opacity: isDragging ? 0.5 : 1 }}
        />
    );
}

const SpatialLayerDialogComponent = observer(() => {
    const chart = useChart<SpatialDataMdvReactConfig, SpatialDataMdvReact>();
    const config = useConfig<SpatialDataMdvReactConfig>();
    const { spatialData } = useChartSpatialData();
    const stack = config.renderStack;
    const mutateStack = useRenderStackMutation();
    const [selectedEntryId, setSelectedEntryId] = useState<string | null>(null);
    const [spatialCollapsed, setSpatialCollapsed] = useState(false);
    // Layer groups whose children are hidden in the list. Dialog state only, not saved.
    const [collapsedGroups, setCollapsedGroups] = useState<ReadonlySet<string>>(() => new Set());
    // The row being dragged, the row it is over, and how many levels it has been dragged sideways.
    const [drag, setDrag] = useState<{ activeId: string; overId: string; depthOffset: number } | null>(null);
    // The open add-layer menu: where it hangs, and the group it adds into (null: the top level).
    const [addMenu, setAddMenu] = useState<{ anchor: HTMLElement; parentId: string | null } | null>(null);
    const [deckCollapsed, setDeckCollapsed] = useState(false);
    const coordinateSystem = //todo datasources schema including spatial regions when referring to dataStore here.
        (chart.dataStore.regions?.all_regions[config.region] as { spatial?: { coordinate_system?: string } })?.spatial
            ?.coordinate_system ?? null;

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
        const options: InsertOption[] = [GROUP_INSERT_OPTION];
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
        const [group, ...elements] = options;
        return [group, ...elements.sort((a, b) => a.label.localeCompare(b.label))];
    }, [coordinateSystem, spatialData, stack]);

    // Memoized by hand: the React Compiler does not optimize observer()-wrapped
    // components — it doesn't recognize mobx-react-lite's `observer` HOC, so every
    // component in this file is invisible to it (verified: a plain sibling panel
    // emits compiler memo-cache output, this all-observer file emits none, in both
    // dev and build). Without explicit useCallback these reallocate on every
    // observer re-render (e.g. per contrast-brush tick), defeating the memoized rows.
    const toggleGroup = useCallback((groupId: string) => {
        setCollapsedGroups((current) => {
            const next = new Set(current);
            if (!next.delete(groupId)) next.add(groupId);
            return next;
        });
    }, []);

    const onDragStart = useCallback((event: DragStartEvent) => {
        const id = String(event.active.id);
        setDrag({ activeId: id, overId: id, depthOffset: 0 });
    }, []);
    const onDragMove = useCallback((event: DragMoveEvent) => {
        const depthOffset = Math.round(event.delta.x / TREE_INDENT_PX);
        setDrag((current) => (current && current.depthOffset !== depthOffset ? { ...current, depthOffset } : current));
    }, []);
    const onDragOver = useCallback((event: DragOverEvent) => {
        const overId = event.over ? String(event.over.id) : null;
        setDrag((current) => (current && overId && current.overId !== overId ? { ...current, overId } : current));
    }, []);
    const onDragCancel = useCallback(() => setDrag(null), []);

    const onUngroup = useCallback(
        (groupId: string) => {
            mutateStack((current) => {
                ungroupRenderStackEntry(current, groupId);
            });
            // the group is gone; its layers keep any failures of their own
            chart.clearLayerFailure(groupId);
            if (selectedEntryId === groupId) setSelectedEntryId(null);
        },
        [chart, mutateStack, selectedEntryId],
    );

    const onInsert = useCallback(
        (option: InsertOption, parentId: string | null) => {
            let id: string | null = null;
            mutateStack((current) => {
                id =
                    option.kind === "group"
                        ? insertRenderStackGroup(current, parentId)
                        : insertSpatialRenderStackEntry(
                              current,
                              option.type,
                              option.elementKey,
                              defaultPropsForSpatialElement(option.type),
                              parentId,
                          );
            });
            if (parentId) {
                // show what was just added to it
                setCollapsedGroups((current) => {
                    if (!current.has(parentId)) return current;
                    const next = new Set(current);
                    next.delete(parentId);
                    return next;
                });
            }
            if (id) setSelectedEntryId(id);
        },
        [mutateStack],
    );
    const onAddInto = useCallback((groupId: string, anchor: HTMLElement) => {
        setAddMenu({ anchor, parentId: groupId });
    }, []);

    const onRemove = useCallback(
        (id: string) => {
            mutateStack((current) => removeRenderStackEntry(current, id));
            chart.clearLayerFailure(id);
            if (selectedEntryId === id) setSelectedEntryId(null);
        },
        [chart, mutateStack, selectedEntryId],
    );

    if (!stack?.entries.length) {
        return (
            <div className="block w-full p-3">
                <Typography color="text.secondary">Layer stack is loading…</Typography>
            </div>
        );
    }

    // Top-level entries other than MDV overlays, and groups other than ones holding only MDV
    // overlays (the data points' group): those are listed with the overlays they hold.
    const spatialEntries = stack.entries.filter((entry) =>
        entry.kind === "group"
            ? !Array.from(walkRenderStackEntries(entry.children)).some((child) => child.kind === "host")
            : !isDeckHostLayerId(entry.id),
    );
    // While a group is dragged its children go with it, so they are not rows of their own.
    const treeRows = flattenLayerTree(spatialEntries, collapsedGroups, drag?.activeId ?? null);
    const projected = drag && projectLayerTreeDrop(treeRows, drag.activeId, drag.overId, drag.depthOffset);
    const onDragEnd = (event: DragEndEvent) => {
        const state = drag;
        setDrag(null);
        if (!state || !event.over) return;
        // from the release itself: the last move may not have reached `drag` yet
        const depthOffset = Math.round(event.delta.x / TREE_INDENT_PX);
        const drop = projectLayerTreeDrop(treeRows, state.activeId, String(event.over.id), depthOffset);
        if (!drop) return;
        const from = findRenderStackSiblings(stack, state.activeId);
        const unchanged =
            from &&
            (from.parent?.id ?? null) === drop.parentId &&
            (from.siblings[from.index - 1]?.id ?? null) === drop.afterId;
        if (unchanged) return;
        mutateStack((current) => {
            moveRenderStackEntry(current, state.activeId, drop.parentId, drop.afterId);
        });
    };
    // MDV overlays always draw above SpatialData layers, so they are listed apart, in a fixed
    // order, and cannot be dragged among them. They may sit inside a group.
    const hasEntry = (deckId: DeckOverlayId) => findRenderStackEntry(stack, deckHostLayerId(deckId)) !== undefined;
    const listedOverlays = LISTED_DECK_OVERLAYS.filter((row) =>
        row.kind === "overlay" ? hasEntry(row.deckId) : row.children.some(hasEntry),
    );

    return (
        <div className="flex flex-col w-full h-full">
            {/* Layer tree */}
            <div
                className="flex-none overflow-y-auto border-b border-gray-200 dark:border-gray-700"
                style={{ maxHeight: "40%" }}
            >
                <>
                    <div className="flex items-center pr-1">
                        <button
                            type="button"
                            className="flex grow items-center gap-1 px-1 py-0.5 text-left hover:bg-gray-100 dark:hover:bg-gray-800"
                            onClick={() => setSpatialCollapsed((c) => !c)}
                        >
                            {spatialCollapsed ? (
                                <ChevronRightIcon fontSize="small" />
                            ) : (
                                <ExpandMoreIcon fontSize="small" />
                            )}
                            <Typography
                                variant="caption"
                                color="text.secondary"
                                className="font-semibold uppercase tracking-wide select-none"
                            >
                                Spatial elements
                            </Typography>
                        </button>
                        <IconButton
                            size="small"
                            onClick={(e) => setAddMenu({ anchor: e.currentTarget, parentId: null })}
                            aria-label="add a layer"
                            title="Add a layer or group"
                        >
                            <AddIcon fontSize="small" />
                        </IconButton>
                    </div>
                    <Collapse in={!spatialCollapsed}>
                        {spatialEntries.length === 0 && (
                            <Typography variant="body2" color="text.secondary" className="px-3 py-1">
                                No SpatialData layers. Add one with +.
                            </Typography>
                        )}
                        <DndContext
                            sensors={sensors}
                            collisionDetection={closestCenter}
                            onDragStart={onDragStart}
                            onDragMove={onDragMove}
                            onDragOver={onDragOver}
                            onDragEnd={onDragEnd}
                            onDragCancel={onDragCancel}
                        >
                            <SortableContext
                                items={treeRows.map((row) => row.id)}
                                strategy={verticalListSortingStrategy}
                            >
                                {treeRows.map((row) => (
                                    <SortableLayerTreeRow
                                        key={row.id}
                                        entryId={row.id}
                                        depth={row.id === drag?.activeId && projected ? projected.depth : row.depth}
                                        selected={selectedEntryId === row.id}
                                        onSelect={setSelectedEntryId}
                                        onRemove={onRemove}
                                        onUngroup={onUngroup}
                                        onAddInto={onAddInto}
                                        expander={
                                            row.isGroup
                                                ? {
                                                      expanded: !collapsedGroups.has(row.id),
                                                      onToggle: () => toggleGroup(row.id),
                                                  }
                                                : undefined
                                        }
                                    />
                                ))}
                            </SortableContext>
                        </DndContext>
                    </Collapse>
                </>
                {listedOverlays.length > 0 && (
                    <>
                        <button
                            type="button"
                            className="flex w-full items-center gap-1 px-1 py-0.5 text-left hover:bg-gray-100 dark:hover:bg-gray-800"
                            onClick={() => setDeckCollapsed((c) => !c)}
                        >
                            {deckCollapsed ? (
                                <ChevronRightIcon fontSize="small" />
                            ) : (
                                <ExpandMoreIcon fontSize="small" />
                            )}
                            <Typography
                                variant="caption"
                                color="text.secondary"
                                className="font-semibold uppercase tracking-wide select-none"
                            >
                                MDV data
                            </Typography>
                        </button>
                        <Collapse in={!deckCollapsed}>
                            {listedOverlays.map((row) =>
                                row.kind === "overlay" ? (
                                    <LayerTreeRow
                                        key={row.deckId}
                                        entryId={deckHostLayerId(row.deckId)}
                                        selected={selectedEntryId === deckHostLayerId(row.deckId)}
                                        onSelect={setSelectedEntryId}
                                        onRemove={onRemove}
                                    />
                                ) : (
                                    <OverlayGroupRow
                                        key={row.label}
                                        label={row.label}
                                        deckIds={row.children}
                                        selectedEntryId={selectedEntryId}
                                        onSelect={setSelectedEntryId}
                                        onRemove={onRemove}
                                    />
                                ),
                            )}
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
                        <LayerHeader entryId={selectedEntryId} />
                        <LayerFailureDetails entryId={selectedEntryId} />
                        <LayerStatusDetails entryId={selectedEntryId} />
                        <LayerDetails entryId={selectedEntryId} />
                        <LayerBlendModePanel entryId={selectedEntryId} />
                        <GroupLayerPanel entryId={selectedEntryId} onUngroup={onUngroup} />
                    </ErrorBoundary>
                ) : (
                    <Typography color="text.secondary" variant="body2">
                        Select a layer to view its properties.
                    </Typography>
                )}
            </div>

            <Menu
                anchorEl={addMenu?.anchor}
                open={addMenu !== null}
                onClose={() => setAddMenu(null)}
                slotProps={{ paper: { sx: { maxHeight: 360 } } }}
            >
                {insertOptions.flatMap((option, index) => {
                    const item = (
                        <MenuItem
                            key={option.label}
                            dense
                            onClick={() => {
                                const parentId = addMenu?.parentId ?? null;
                                setAddMenu(null);
                                onInsert(option, parentId);
                            }}
                        >
                            {option.kind === "group" ? "Group (drag layers into it)" : option.label}
                        </MenuItem>
                    );
                    return index === 1
                        ? [
                              <Divider key="divider" />,
                              <ListSubheader key="elements">SpatialData elements</ListSubheader>,
                              item,
                          ]
                        : [item];
                })}
            </Menu>
        </div>
    );
});

export default SpatialLayerDialogComponent;
