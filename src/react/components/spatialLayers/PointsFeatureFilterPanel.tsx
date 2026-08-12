import { Alert, Button, Checkbox, FormControlLabel, TextField, Typography } from "@mui/material";
import { featureNamesForCodes, isPointsWorkerEnabled, resolveFeatureSelectionCodes } from "@spatialdata/core";
import { featureCodeToRgb } from "@spatialdata/layers";
import { describeFeatureRowState, featureRowOpacity, usePointsFeatureState } from "@spatialdata/vis";
import { memo, useEffect, useMemo, useRef, useState } from "react";

import {
    type PointsFeatureFilterConfig,
    type PointsLayerUpdate,
    samePointsFeatureFilterConfig,
} from "@/react/spatialdata/points_layer_config";

/** Above this many features the list gets a search box; below it, scrolling is enough. */
const FEATURE_LIST_SEARCH_THRESHOLD = 100;

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
    const effectiveRgb = (name: string, code: number): [number, number, number] =>
        colorOverrides?.[name] ?? featureCodeToRgb(code);

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
            if (pending) updateLayer({ featureColorOverrides: pending });
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
    // necessary but not sufficient: the scan that uses it runs in the core points
    // worker, and `loadPointsMatchingFeatureCodes` throws outright without one rather
    // than falling back to the main thread. `ensurePointsWorker` starts it, so this is
    // normally true — but if that ever fails the row must not invite a click that
    // cannot work, which is the state this app shipped in before core 0.8.0.
    const canScanOnDemand = supportsOnDemandLoad && isPointsWorkerEnabled();
    /** The element could scan, but the worker it needs never started. */
    const workerBlocksScan = supportsOnDemandLoad && !canScanOnDemand;
    // ONE classification pass over the catalog, not one per consumer. This was a
    // `rowInfo(code)` function called from two `reduce`s and again from the row map, so
    // `describeFeatureRowState` ran three times per feature — a ~36ms floor on every
    // render at 541 features, paid even when the search box narrowed the list to one row.
    // Deliberately not `useMemo`d: `loadedMatchingCodes` is a fresh Set per engine read,
    // so a dependency array would miss every time and only add the cost of checking.
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

            {notLoadedCount > 0 ? (
                <Alert severity={canScanOnDemand ? "info" : "warning"} sx={{ py: 0 }}>
                    <Typography variant="caption">
                        {notLoadedCount} of {entries.length} feature
                        {entries.length === 1 ? "" : "s"}{" "}
                        {canScanOnDemand
                            ? "not loaded yet (greyed below) — selecting one loads it on demand."
                            : workerBlocksScan
                              ? "not in the loaded sample (greyed below). The points worker that fetches them did not start, so raise the memory cap to bring more in."
                              : "not in the loaded sample (greyed below). This dataset has no feature index, so they can't be shown until the memory cap is raised or it's rewritten with one."}
                    </Typography>
                </Alert>
            ) : null}

            {partialCount > 0 ? (
                <Typography variant="caption" color="text.secondary">
                    {partialCount} of {entries.length} feature{entries.length === 1 ? "" : "s"} only partly loaded — the
                    resident window is capped, so the canvas is drawing a sample of each. Select one to fetch it in
                    full, or raise the memory cap.
                </Typography>
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

            <div className="grid max-h-56 gap-0.5 overflow-y-auto py-1">
                {visibleEntries.map((entry) => {
                    const { resident, rendered, selected, state } = rowInfo(entry.code);
                    const overridden = colorOverrides?.[entry.name] !== undefined;
                    const rgb = effectiveRgb(entry.name, entry.code);
                    // Once dataset totals land, keep showing the resident tally too when
                    // it falls short. Dropping it is what let a capped element print
                    // "1,182,402" beside a row drawing a tenth of that. `partial` already
                    // excludes features a scan has since supplied whole.
                    const shortfall = state.tone === "partial" ? residentShortfall(entry) : undefined;
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
                            key={entry.code}
                            sx={{ ...rowLabelSx, opacity: featureRowOpacity(state) }}
                            title={title}
                            onMouseEnter={() => setHighlightedFeature(entry.code)}
                            onMouseLeave={() => setHighlightedFeature(null)}
                            control={
                                <Checkbox
                                    size="small"
                                    sx={checkboxSx}
                                    checked={selected}
                                    onChange={(event) => toggleFeature(entry.code, event.target.checked)}
                                />
                            }
                            label={
                                <span className="flex items-center gap-1.5">
                                    <FeatureColorSwatch
                                        name={entry.name}
                                        rgb={rgb}
                                        overridden={overridden}
                                        onPick={(next) => setColorOverride(entry.name, next)}
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
                                                clearColorOverride(entry.name);
                                            }}
                                        >
                                            ⟲
                                        </button>
                                    ) : null}
                                    {hasAnyCounts ? (
                                        <Typography
                                            variant="caption"
                                            color="text.secondary"
                                            className="ml-auto shrink-0"
                                            title={
                                                shortfall !== undefined
                                                    ? `${shortfall.toLocaleString()} of ${entry.count?.toLocaleString()} points are inside the memory cap`
                                                    : countIsPartial(entry)
                                                      ? "Points loaded so far (resident window) — dataset total still counting"
                                                      : "Points in the dataset"
                                            }
                                        >
                                            {shortfall !== undefined ? (
                                                <>
                                                    {/* The resident figure is what is on
                                                        screen, so it carries the emphasis;
                                                        the dataset total is context. */}
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
                                                    {countIsPartial(entry) ? "≥" : ""}
                                                    {effectiveCount(entry)?.toLocaleString() ?? "—"}
                                                </>
                                            )}
                                        </Typography>
                                    ) : null}
                                </span>
                            }
                        />
                    );
                })}
                {showSearch && visibleEntries.length === 0 ? (
                    <Typography variant="caption" color="text.secondary">
                        No features match your search.
                    </Typography>
                ) : null}
            </div>
        </div>
    );
}

/**
 * Memoised on the three config fields it reads, so the cosmetic edits that wake the
 * layer dialog — opacity, point size, memory cap — do not re-render 541 feature rows.
 * The comparator lives with the type it compares; see `PointsFeatureFilterConfig`.
 *
 * Engine-driven updates are unaffected: `usePointsFeatureState` subscribes this
 * component through `useSyncExternalStore`, which re-renders it regardless of props.
 */
export default memo(PointsFeatureFilterPanel, (prev, next) => {
    return prev.updateLayer === next.updateLayer && samePointsFeatureFilterConfig(prev.config, next.config);
});
