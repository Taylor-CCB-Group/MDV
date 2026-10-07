import { Alert, Button, Checkbox, FormControlLabel, TextField, Typography } from "@mui/material";
import { featureNamesForCodes, isParquetWorkerEnabled, resolveFeatureSelectionCodes } from "@spatialdata/core";
import { featureCodeToRgb } from "@spatialdata/layers";
import { describeFeatureRowState, featureRowOpacity, usePointsFeatureState } from "@spatialdata/vis";
import { useVirtualizer } from "@tanstack/react-virtual";
import { type ReactNode, memo, useEffect, useMemo, useRef, useState } from "react";

import {
    type PointsFeatureFilterConfig,
    type PointsLayerUpdate,
    samePointsFeatureFilterConfig,
} from "@/react/spatialdata/points_layer_config";

/** Above this many features the list gets a search box; below it, scrolling is enough. */
const FEATURE_LIST_SEARCH_THRESHOLD = 100;

/**
 * Starting guess for a row's height: ~22px for a caption-sized row carrying a small MUI
 * checkbox, plus the 2px that used to be the list grid's `gap-0.5`. Only the length of
 * the scrollbar before a row has been seen depends on it — every mounted row reports its
 * real height back through `measureElement`, so a wrong guess settles itself as you
 * scroll instead of overlapping rows, which a hard-coded row height would.
 */
const FEATURE_ROW_HEIGHT = 24;

/** Rows kept mounted either side of the window, so a flick-scroll doesn't show blanks. */
const FEATURE_ROW_OVERSCAN = 6;

/**
 * `FormControlLabel` rather than a bare `<label>` around a MUI `Checkbox`: the
 * checkbox's real `<input>` is nested inside a component, which an a11y linter can't
 * see, and this is the pairing MUI intends anyway. Its defaults are built for form
 * rows, so strip the margin and let the label take the remaining width.
 */
const rowLabelSx = {
    m: 0,
    gap: 0.75,
    minWidth: 0,
    "& .MuiFormControlLabel-label": { flex: 1, minWidth: 0 },
} as const;

const checkboxSx = { p: 0.25 } as const;

/** Dotted underline so the terse coverage counts read as "there is more here". */
const HINT_CLASS = "cursor-help underline decoration-dotted underline-offset-2";

const NOT_LOADED_HINT =
    "Greyed below: no points inside the memory cap. Selecting one fetches it on demand, or raise the cap under Advanced.";

const PARTIAL_HINT =
    "Drawn from a sample — some of their points are outside the memory cap. Select one to fetch it in full, or raise the cap under Advanced.";

const pluralFeatures = (count: number) => `${count} feature${count === 1 ? "" : "s"}`;

type Props = {
    config: PointsFeatureFilterConfig;
    updateLayer: PointsLayerUpdate;
};

type CatalogEntry = { code: number; name: string; count?: number };

type RowState = {
    resident: boolean;
    rendered: boolean;
    selected: boolean;
    state: ReturnType<typeof describeFeatureRowState>;
};

/** For a code the catalog does not list — unreachable from the rendered rows, which
 * are drawn from `entries`, but a lookup must return something. */
const UNKNOWN_ROW_STATE: RowState = {
    resident: false,
    rendered: false,
    selected: false,
    state: describeFeatureRowState({
        resident: false,
        rendered: false,
        selected: false,
        scanning: false,
        supportsOnDemandLoad: false,
        residentKnown: false,
    }),
};

const hex2 = (value: number) => Math.max(0, Math.min(255, value)).toString(16).padStart(2, "0");

const rgbToHex = ([r, g, b]: readonly [number, number, number]) => `#${hex2(r)}${hex2(g)}${hex2(b)}`;

const hexToRgb = (hex: string): [number, number, number] => {
    const n = Number.parseInt(hex.slice(1), 16);
    return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
};

/**
 * The swatch IS the picker: the span shows the feature's effective colour and a
 * transparent native colour input sits over it. `inline-block` + `border-box` keep
 * the 12×12 box (and its border) exact inside the flex row.
 */
function FeatureColorSwatch({
    name,
    rgb,
    overridden,
    onPick,
}: {
    name: string;
    rgb: [number, number, number];
    overridden: boolean;
    onPick: (rgb: [number, number, number]) => void;
}) {
    return (
        <span
            className="relative inline-block box-border h-3 w-3 shrink-0 rounded-sm"
            style={{
                background: `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`,
                border: overridden ? "1px solid #6cb6ff" : "1px solid rgba(128, 128, 128, 0.5)",
                boxShadow: overridden ? "0 0 0 1px #6cb6ff" : undefined,
            }}
            title={`${name} colour${overridden ? " (overridden)" : ""}`}
        >
            <input
                type="color"
                aria-label={`${name} colour`}
                value={rgbToHex(rgb)}
                className="absolute inset-0 m-0 h-full w-full cursor-pointer appearance-none border-none p-0 opacity-0"
                onClick={(event) => event.stopPropagation()}
                onChange={(event) => onPick(hexToRgb(event.target.value))}
            />
        </span>
    );
}

/**
 * One feature row. `memo` is the point of it: `FeatureRowList` re-renders on every scroll
 * event, not just when the window moves, and MUI's `FormControlLabel` + `Checkbox` cost
 * enough per row that re-rendering all ~22 of them each time was ~45ms — still visible as
 * lag on a fast scroll. Every prop below is a primitive, an entry from the catalog, or a
 * value the panel already holds by stable reference, so a scroll that does not change a
 * row skips it entirely and only rows entering the window actually render.
 *
 * That is also the contract to keep: pass derived values, never freshly-built objects or
 * per-row closures — either turns the shallow compare into a guaranteed miss and quietly
 * puts the ~45ms back.
 */
const FeatureRow = memo(function FeatureRow({
    entry,
    row,
    colorOverride,
    shortfall,
    displayCount,
    countIsPartial,
    showCount,
    scanning,
    workerBlocksScan,
    onToggle,
    onHighlight,
    onPickColor,
    onClearColor,
}: {
    entry: CatalogEntry;
    row: RowState;
    /** This feature's colour override, if it has one. Also what makes the row "overridden". */
    colorOverride?: [number, number, number];
    /** Resident points, when the row is drawn from fewer points than the dataset holds. */
    shortfall?: number;
    /** The count to print: the dataset total, or the running resident tally before it lands. */
    displayCount?: number;
    /** `displayCount` is still only the running tally, so it prints with a "≥". */
    countIsPartial: boolean;
    /** No counts anywhere in the catalog yet, so no count column on any row. */
    showCount: boolean;
    scanning: boolean;
    workerBlocksScan: boolean;
    onToggle: (code: number, checked: boolean) => void;
    onHighlight: (code: number | null) => void;
    onPickColor: (name: string, rgb: [number, number, number]) => void;
    onClearColor: (name: string) => void;
}) {
    const { resident, rendered, selected, state } = row;
    const overridden = colorOverride !== undefined;
    const rgb = colorOverride ?? featureCodeToRgb(entry.code);
    const countStr = entry.count !== undefined ? ` · ${entry.count.toLocaleString()} pts` : "";
    // Multi-line diagnostic: the human state and its reason, then the raw
    // signals that drove it (what made this row grey, or not).
    const title = `${entry.name} · code ${entry.code}${countStr}\n${state.label}: ${state.reason}\n${
        workerBlocksScan
            ? "(This element does have a feature index; the points worker it needs can't start — SpatialData.js#148.)\n"
            : ""
    }[resident=${resident ? "y" : "n"} rendered=${rendered ? "y" : "n"} selected=${selected ? "y" : "n"} scan=${scanning ? "running" : "idle"}]`;
    return (
        <FormControlLabel
            // `width` because the old grid container stretched its items and
            // `FeatureRowList`'s row wrapper does not; without it the label
            // shrinks to its text and the right-aligned point count collapses
            // onto the name.
            sx={{ ...rowLabelSx, width: "100%", opacity: featureRowOpacity(state) }}
            title={title}
            onMouseEnter={() => onHighlight(entry.code)}
            onMouseLeave={() => onHighlight(null)}
            // Focus too, so tabbing through the list highlights on the
            // canvas the same way hovering does. These bubble from the
            // row's checkbox, which is the focusable element.
            onFocus={() => onHighlight(entry.code)}
            onBlur={() => onHighlight(null)}
            control={
                <Checkbox
                    size="small"
                    sx={checkboxSx}
                    checked={selected}
                    onChange={(event) => onToggle(entry.code, event.target.checked)}
                />
            }
            label={
                <span className="flex items-center gap-1.5">
                    <FeatureColorSwatch
                        name={entry.name}
                        rgb={rgb}
                        overridden={overridden}
                        onPick={(next) => onPickColor(entry.name, next)}
                    />
                    <Typography variant="caption" noWrap>
                        {entry.name}
                    </Typography>
                    {overridden ? (
                        <button
                            type="button"
                            title="Reset to default colour"
                            className="shrink-0 cursor-pointer rounded border border-[hsl(var(--border))] px-1 text-[11px] leading-tight text-[hsl(var(--muted-foreground))]"
                            onClick={(event) => {
                                event.stopPropagation();
                                event.preventDefault();
                                onClearColor(entry.name);
                            }}
                        >
                            ⟲
                        </button>
                    ) : null}
                    {showCount ? (
                        <Typography
                            variant="caption"
                            color="text.secondary"
                            className="ml-auto shrink-0"
                            title={
                                shortfall !== undefined
                                    ? `${shortfall.toLocaleString()} of ${entry.count?.toLocaleString()} points are inside the memory cap`
                                    : countIsPartial
                                      ? "Points loaded so far (resident window) — dataset total still counting"
                                      : "Points in the dataset"
                            }
                        >
                            {shortfall !== undefined ? (
                                <>
                                    {/* The resident figure is what is on screen, so it carries
                                        the emphasis; the dataset total is context. */}
                                    <span className="text-[hsl(var(--warning,38_92%_50%))]">
                                        {shortfall.toLocaleString()}
                                    </span>
                                    <span className="opacity-60">
                                        {" / "}
                                        {entry.count?.toLocaleString()}
                                    </span>
                                </>
                            ) : (
                                <>
                                    {countIsPartial ? "≥" : ""}
                                    {displayCount?.toLocaleString() ?? "—"}
                                </>
                            )}
                        </Typography>
                    ) : null}
                </span>
            }
        />
    );
});

/**
 * The virtualised feature list — mount only the rows the fixed-height scroll box can
 * actually show. A Xenium transcripts element has ~12,400 features, and one
 * `FormControlLabel` + `Checkbox` each put ~88,000 nodes on the page and blocked the main
 * thread for long enough that the tab stopped answering at all: the panel was unusable on
 * exactly the datasets the feature filter exists for.
 *
 * Its own component, and NOT inlinable back into the panel: `useVirtualizer` re-renders
 * whoever calls it on every scroll frame, and the panel's render is O(catalog) — it
 * classifies all ~12,400 features to count how many are unloaded. Called from there, each
 * scroll frame cost ~310ms and rows arrived visibly after the scrollbar had moved. Here
 * the hook re-renders the window's rows and nothing else; `renderRow` closes over whatever
 * the panel's last render computed, so scrolling never re-derives it.
 *
 * `lint:react-compiler` warns "incompatible library" on `useVirtualizer` — it returns
 * functions the compiler cannot memoise safely, so it skips this component. That is the
 * outcome we want anyway: memoising a list whose rows depend on mutable engine state is
 * how it would go stale.
 */
function FeatureRowList({
    count,
    getItemKey,
    renderRow,
    children,
}: {
    count: number;
    /** The row's feature code — the identity `measureElement` caches its height under. */
    getItemKey: (index: number) => number;
    renderRow: (index: number) => ReactNode;
    /** Drawn after the rows, inside the scroll box: the "nothing matched" message. */
    children?: ReactNode;
}) {
    const scrollRef = useRef<HTMLDivElement>(null);
    const rowVirtualizer = useVirtualizer({
        count,
        getScrollElement: () => scrollRef.current,
        estimateSize: () => FEATURE_ROW_HEIGHT,
        overscan: FEATURE_ROW_OVERSCAN,
        // Key measurements by feature, not by position: the search box reshuffles which
        // entry sits at which index, and an index-keyed cache would hand a row the
        // measurement of whichever feature happened to be there before.
        getItemKey,
    });

    return (
        <div ref={scrollRef} className="max-h-56 overflow-y-auto py-1">
            {/* Spacer as tall as the whole list, with only the scrolled-to window of rows
                inside it. Absolute positioning has no flow for the old grid's `gap-0.5` to
                act on, so the gap moves into each row's own padding — where it stays part
                of what `measureElement` measures. */}
            <div className="relative w-full" style={{ height: rowVirtualizer.getTotalSize() }}>
                {rowVirtualizer.getVirtualItems().map((virtualItem) => (
                    <div
                        key={virtualItem.key}
                        // `data-index` is how the default `measureElement` finds which row it
                        // just measured; without it every row measures as index 0.
                        data-index={virtualItem.index}
                        ref={rowVirtualizer.measureElement}
                        className="absolute left-0 top-0 w-full pb-0.5"
                        style={{ transform: `translateY(${virtualItem.start}px)` }}
                    >
                        {renderRow(virtualItem.index)}
                    </div>
                ))}
            </div>
            {children}
        </div>
    );
}

function PointsFeatureFilterPanel({ config, updateLayer }: Props) {
    // Opt out of the React Compiler. `usePointsFeatureState` re-renders this component
    // on every engine notify (via useSyncExternalStore) but reads mutable engine state
    // the compiler cannot see as a dependency, so it would memoize this JSX and hold the
    // pre-catalog branch on screen after the catalog arrives. Scoped to this leaf.
    "use no memo";
    const {
        catalog,
        catalogLoading,
        catalogRefining,
        residentCodes,
        loadedMatchingCodes,
        supportsOnDemandLoad,
        matchingLoadState,
        residentFeatureCounts,
        requestCatalog,
        setHighlightedFeature,
        retryFailedLoads,
    } = usePointsFeatureState(config);

    const [searchQuery, setSearchQuery] = useState("");

    // Ask for the full-dataset catalog whenever the panel is shown for a layer. The
    // engine dedupes, so this only upgrades the instant resident-subset preview to the
    // complete list plus counts.
    useEffect(() => {
        requestCatalog();
    }, [requestCatalog]);

    // Drop any lingering hover emphasis when the panel unmounts or its layer changes,
    // so a highlight doesn't stick on the canvas after the pointer is long gone.
    useEffect(() => () => setHighlightedFeature(null), [setHighlightedFeature]);

    const entries = useMemo<CatalogEntry[]>(() => catalog?.entries ?? [], [catalog?.entries]);
    const hasCounts = entries.some((entry) => entry.count !== undefined);
    // Dataset totals only arrive with the catalog's counts scan. Until then fall back to
    // the running tally over the resident window, marked "≥" so a partial figure is never
    // mistaken for a total.
    const partialCounts = residentFeatureCounts;
    const hasAnyCounts = hasCounts || (partialCounts?.size ?? 0) > 0;
    const effectiveCount = (entry: CatalogEntry) => entry.count ?? partialCounts?.get(entry.code);
    const countIsPartial = (entry: CatalogEntry) =>
        entry.count === undefined && partialCounts?.get(entry.code) !== undefined;
    /** Resident points for a feature when that is meaningfully LESS than the dataset —
     * a shortfall worth showing. `undefined` otherwise. */
    const residentShortfall = (entry: CatalogEntry): number | undefined => {
        if (entry.count === undefined) return undefined;
        const resident = residentFeatureCounts?.get(entry.code);
        if (resident === undefined || resident >= entry.count) return undefined;
        return resident;
    };

    // The selection persists as NAMES (see PointsLayerConfig.featureNames); everything
    // below works in codes, resolved once against the catalog already being rendered.
    const selection = resolveFeatureSelectionCodes(config, catalog);
    const allSelected = selection === undefined;
    const noneSelected = selection !== undefined && selection.length === 0;
    const selectedCodes = allSelected ? new Set(entries.map((entry) => entry.code)) : new Set(selection ?? []);

    const sortedEntries = useMemo(() => {
        const list = [...entries];
        const rank = (entry: CatalogEntry) => entry.count ?? partialCounts?.get(entry.code) ?? -1;
        if (hasAnyCounts) {
            list.sort((left, right) => rank(right) - rank(left) || left.name.localeCompare(right.name));
        } else {
            list.sort((left, right) => left.name.localeCompare(right.name));
        }
        return list;
    }, [entries, hasAnyCounts, partialCounts]);

    const visibleEntries = useMemo(() => {
        const query = searchQuery.trim().toLowerCase();
        if (!query) return sortedEntries;
        return sortedEntries.filter((entry) => entry.name.toLowerCase().includes(query));
    }, [sortedEntries, searchQuery]);

    // Write NAMES, and clear any legacy `featureCodes` so the two cannot disagree —
    // names win when both are set, and a stale code list in a saved config is exactly
    // the confusion names exist to remove.
    const setSelectedCodes = (nextCodes: number[] | undefined) => {
        updateLayer({
            featureNames: nextCodes ? featureNamesForCodes(nextCodes, catalog) : undefined,
            featureCodes: undefined,
        });
    };

    const colorOverrides = config.featureColorOverrides;

    // `<input type="color">` fires change continuously while the picker is dragged, and
    // each commit is a config write → new palette → deck layer update on a layer that can
    // hold millions of points. Coalesce to one write per frame: the canvas still previews
    // live, but the work is bounded by the display rather than by event rate.
    // The pending value is the FULL next overrides map, not one entry, so two features
    // recoloured inside a frame both survive, and the merge base is taken at schedule
    // time rather than depending on which render created the handler.
    const pendingColorRef = useRef<Record<string, [number, number, number]> | null>(null);
    const colorFrameRef = useRef<number | null>(null);
    useEffect(
        () => () => {
            if (colorFrameRef.current !== null) cancelAnimationFrame(colorFrameRef.current);
        },
        [],
    );
    const setColorOverride = (name: string, rgb: [number, number, number]) => {
        pendingColorRef.current = {
            ...(pendingColorRef.current ?? colorOverrides ?? {}),
            [name]: rgb,
        };
        if (colorFrameRef.current !== null) return;
        colorFrameRef.current = requestAnimationFrame(() => {
            colorFrameRef.current = null;
            const pending = pendingColorRef.current;
            pendingColorRef.current = null;
            // Clearing the last override inside the same frame as a recolour leaves an
            // empty map here. Writing `{}` would persist an override record that says
            // "no overrides" — `undefined` is how the rest of this file spells that.
            if (pending) {
                updateLayer({
                    featureColorOverrides: Object.keys(pending).length > 0 ? pending : undefined,
                });
            }
        });
    };
    const clearColorOverride = (name: string) => {
        // A coalesced write may still be queued, and it carries the whole map — letting
        // it land after the clear would put the override straight back.
        if (pendingColorRef.current && name in pendingColorRef.current) {
            delete pendingColorRef.current[name];
        }
        if (!colorOverrides || !(name in colorOverrides)) return;
        const next = { ...colorOverrides };
        delete next[name];
        updateLayer({
            featureColorOverrides: Object.keys(next).length > 0 ? next : undefined,
        });
    };

    const toggleFeature = (code: number, checked: boolean) => {
        const current = new Set(allSelected ? entries.map((entry) => entry.code) : (selection ?? []));
        if (checked) current.add(code);
        else current.delete(code);
        if (current.size === 0) return setSelectedCodes([]);
        if (current.size === entries.length) return setSelectedCodes(undefined);
        setSelectedCodes([...current].sort((left, right) => left - right));
    };

    // Only block on loading when there is nothing to show at all. The catalog scan
    // publishes names/codes before its slower per-feature counts pass, so the list is
    // usable — selectable and colourable — while the counts column fills in.
    if (catalogLoading && !catalog) {
        return (
            <Typography variant="caption" color="text.secondary">
                Loading features…
            </Typography>
        );
    }

    if (catalog === undefined) {
        return (
            <div className="grid gap-2 justify-items-start">
                <Typography variant="caption" color="text.secondary">
                    Feature list not loaded.
                </Typography>
                <Button size="small" variant="outlined" onClick={() => requestCatalog()}>
                    Load feature list
                </Button>
            </div>
        );
    }

    if (!catalog || entries.length === 0) {
        return (
            <Typography variant="caption" color="text.secondary">
                {catalog === null
                    ? "No feature catalog for this points layer (no feature_key, or an encoding this dataset size doesn't support)."
                    : "No features found in the feature catalog."}
            </Typography>
        );
    }

    const selectedCount = noneSelected ? 0 : allSelected ? entries.length : selectedCodes.size;
    const showSearch = entries.length > FEATURE_LIST_SEARCH_THRESHOLD;
    // A feature counts as on screen if it is in the instant resident preview OR in the
    // last-completed feature-index scan. Keying off what is rendered — not the current
    // scan's settled state — keeps already-loaded rows un-greyed while a newly added
    // feature's scan is still in flight.
    const residentKnown = residentCodes !== undefined;
    const scanning = matchingLoadState?.loading ?? false;
    // `supportsOnDemandLoad` answers "does this element have a feature index", which is
    // necessary but not sufficient: the scan that uses it runs in the core parquet
    // worker, and `loadPointsMatchingFeatureCodes` throws outright without one rather
    // than falling back to the main thread. `ensureParquetWorker` starts it, so this is
    // normally true — but a worker that fails to load is now detected and switched off,
    // so `isParquetWorkerEnabled()` reports `false` rather than staying optimistic, and
    // the row must not invite a click that cannot work.
    const canScanOnDemand = supportsOnDemandLoad && isParquetWorkerEnabled();
    /** The element could scan, but the worker it needs never started. */
    const workerBlocksScan = supportsOnDemandLoad && !canScanOnDemand;
    // ONE classification pass over the catalog, not one per consumer. This was a
    // `rowInfo(code)` function called from two `reduce`s and again from the row map, so
    // `describeFeatureRowState` ran three times per feature — a ~36ms floor on every
    // render at 541 features, paid even when the search box narrowed the list to one row.
    // Deliberately not `useMemo`d: `loadedMatchingCodes` is a fresh Set per engine read,
    // so a dependency array would miss every time and only add the cost of checking.
    //
    // Virtualising the list did not make this pass optional, only the reason for it
    // narrower. The rows no longer drive it — a dozen-odd of them are mounted — but
    // `notLoadedCount` and `partialCount` are counts over the WHOLE catalog, so every
    // entry still has to be classified once. The window's `rowInfo` lookups now ride
    // along on a pass that has to happen anyway.
    //
    // This is what makes rendering this component expensive (~310ms at 12,448 features),
    // and so what `FeatureRowList` exists to keep off the scroll path.
    const rowStates = new Map<number, RowState>();
    let notLoadedCount = 0;
    // `partialCount` is drawn-but-incomplete — the opposite failure of understanding to
    // `notLoadedCount`: those rows look entirely healthy, un-greyed with a full dataset
    // count beside them, while most of their points are outside the cap.
    let partialCount = 0;
    for (const entry of entries) {
        const code = entry.code;
        const resident = residentKnown && (residentCodes?.has(code) ?? false);
        const rendered = loadedMatchingCodes?.has(code) ?? false;
        const selected = !noneSelected && (allSelected || selectedCodes.has(code));
        const state = describeFeatureRowState({
            resident,
            rendered,
            selected,
            scanning,
            supportsOnDemandLoad: canScanOnDemand,
            residentKnown,
            residentPointCount: residentFeatureCounts?.get(code),
            datasetPointCount: entry.count,
        });
        rowStates.set(code, { resident, rendered, selected, state });
        // Both counters were guarded on `residentKnown` when they were their own passes;
        // without it every feature reads as not-loaded before the codes land.
        if (residentKnown) {
            if (state.greyed) notLoadedCount += 1;
            if (state.tone === "partial") partialCount += 1;
        }
    }
    const rowInfo = (code: number): RowState => rowStates.get(code) ?? UNKNOWN_ROW_STATE;

    return (
        <div className="grid gap-2">
            <div className="flex items-baseline gap-2">
                <Typography variant="caption">Features ({catalog.featureKey})</Typography>
                <Typography variant="caption" color="text.secondary">
                    {selectedCount}/{entries.length} selected
                    {hasAnyCounts ? (hasCounts ? " · by count" : " · by count so far") : ""}
                </Typography>
            </div>

            {catalogRefining ? (
                <Typography variant="caption" color="text.secondary">
                    Loading the full feature list…
                </Typography>
            ) : null}
            {catalogLoading && !catalogRefining ? (
                // The list is already usable; only the per-feature counts are outstanding.
                <Typography variant="caption" color="text.secondary">
                    Counting features…
                </Typography>
            ) : null}

            {notLoadedCount > 0 || partialCount > 0 ? (
                canScanOnDemand ? (
                    // The healthy case: everything missing is one click from being fetched,
                    // so it is a status line rather than an Alert. The explanation moves to
                    // the tooltips — as prose it was three lines about the memory cap
                    // sitting above the list it was describing.
                    <Typography variant="caption" color="text.secondary">
                        {notLoadedCount > 0 ? (
                            <span className={HINT_CLASS} title={NOT_LOADED_HINT}>
                                {notLoadedCount} not loaded
                            </span>
                        ) : null}
                        {notLoadedCount > 0 && partialCount > 0 ? " · " : null}
                        {partialCount > 0 ? (
                            <span className={HINT_CLASS} title={PARTIAL_HINT}>
                                {partialCount} partly loaded
                            </span>
                        ) : null}
                    </Typography>
                ) : (
                    // Nothing the user can click will fix these, so they keep the Alert.
                    <Alert severity="warning" sx={{ py: 0 }}>
                        <Typography variant="caption">
                            {notLoadedCount > 0
                                ? `${pluralFeatures(notLoadedCount)} not in the loaded sample (greyed below)`
                                : `${pluralFeatures(partialCount)} only partly loaded`}
                            {notLoadedCount > 0 && partialCount > 0
                                ? `, and ${partialCount} more only partly loaded`
                                : ""}
                            .{" "}
                            {workerBlocksScan
                                ? "The points worker that fetches the rest did not start, so raise the memory cap under Advanced to bring more in."
                                : "This dataset has no feature index, so the rest can't be fetched on demand — raise the memory cap under Advanced, or rewrite the element with one."}
                        </Typography>
                    </Alert>
                )
            ) : null}

            {matchingLoadState?.failed ? (
                // A failed scan still DRAWS: the render path falls back to filtering the
                // resident batch, so the canvas shows whichever part of the selection was
                // inside the cap. Saying so matters more than the error text — without it
                // the partial view reads as the complete answer.
                <Alert
                    severity="error"
                    sx={{ py: 0 }}
                    action={
                        matchingLoadState.error?.retryable === false ? undefined : (
                            <Button size="small" onClick={() => retryFailedLoads()}>
                                Retry
                            </Button>
                        )
                    }
                >
                    <Typography variant="caption">
                        Could not load the selected features
                        {matchingLoadState.error ? `: ${matchingLoadState.error.message}` : "."} Showing only the
                        selected points already in memory.
                    </Typography>
                </Alert>
            ) : matchingLoadState ? (
                <Typography variant="caption" color={matchingLoadState.loading ? "primary" : "text.secondary"}>
                    {matchingLoadState.loading
                        ? `Loading selected features… ${matchingLoadState.matchedRows.toLocaleString()} points so far`
                        : matchingLoadState.covered
                          ? `Selection served from ${matchingLoadState.matchedRows.toLocaleString()} points in memory (no re-scan)`
                          : `${matchingLoadState.matchedRows.toLocaleString()} points loaded for this selection`}
                </Typography>
            ) : null}

            <div className="flex items-center gap-4">
                <FormControlLabel
                    sx={rowLabelSx}
                    control={
                        <Checkbox
                            size="small"
                            sx={checkboxSx}
                            checked={allSelected}
                            onChange={(event) => {
                                if (event.target.checked) setSelectedCodes(undefined);
                            }}
                        />
                    }
                    label={<Typography variant="caption">All</Typography>}
                />
                <FormControlLabel
                    sx={rowLabelSx}
                    control={
                        <Checkbox
                            size="small"
                            sx={checkboxSx}
                            checked={noneSelected}
                            onChange={(event) => {
                                if (event.target.checked) setSelectedCodes([]);
                            }}
                        />
                    }
                    label={<Typography variant="caption">None</Typography>}
                />
            </div>

            {showSearch ? (
                <TextField
                    size="small"
                    type="search"
                    placeholder="Search features…"
                    value={searchQuery}
                    onChange={(event) => setSearchQuery(event.target.value)}
                />
            ) : null}

            <FeatureRowList
                count={visibleEntries.length}
                getItemKey={(index) => visibleEntries[index].code}
                renderRow={(index) => {
                    const entry = visibleEntries[index];
                    const row = rowInfo(entry.code);
                    return (
                        <FeatureRow
                            entry={entry}
                            row={row}
                            colorOverride={colorOverrides?.[entry.name]}
                            // Once dataset totals land, keep showing the resident tally too
                            // when it falls short. Dropping it is what let a capped element
                            // print "1,182,402" beside a row drawing a tenth of that.
                            // `partial` already excludes features a scan has since supplied
                            // whole.
                            shortfall={row.state.tone === "partial" ? residentShortfall(entry) : undefined}
                            displayCount={effectiveCount(entry)}
                            countIsPartial={countIsPartial(entry)}
                            showCount={hasAnyCounts}
                            scanning={scanning}
                            workerBlocksScan={workerBlocksScan}
                            onToggle={toggleFeature}
                            onHighlight={setHighlightedFeature}
                            onPickColor={setColorOverride}
                            onClearColor={clearColorOverride}
                        />
                    );
                }}
            >
                {showSearch && visibleEntries.length === 0 ? (
                    <Typography variant="caption" color="text.secondary">
                        No features match your search.
                    </Typography>
                ) : null}
            </FeatureRowList>
        </div>
    );
}

/**
 * Memoised on the three config fields it reads, so the cosmetic edits that wake the
 * layer dialog — opacity, point size, memory cap — do not re-render the panel at all.
 * Virtualisation has since capped the row count, but a re-render still reclassifies
 * every catalog entry for the two coverage counts, so this stays worth having.
 * The comparator lives with the type it compares; see `PointsFeatureFilterConfig`.
 *
 * Engine-driven updates are unaffected: `usePointsFeatureState` subscribes this
 * component through `useSyncExternalStore`, which re-renders it regardless of props.
 */
export default memo(PointsFeatureFilterPanel, (prev, next) => {
    return prev.updateLayer === next.updateLayer && samePointsFeatureFilterConfig(prev.config, next.config);
});
