import { Alert, Button, Checkbox, FormControlLabel, TextField, Typography } from "@mui/material";
import { featureNamesForCodes, isPointsWorkerEnabled, resolveFeatureSelectionCodes } from "@spatialdata/core";
import { featureCodeToRgb } from "@spatialdata/layers";
import { usePointsFeatureState } from "@spatialdata/vis";
import { useEffect, useMemo, useRef, useState } from "react";

import { describeFeatureRowState, featureRowOpacity } from "@/react/spatialdata/points_feature_row_state";
import type { PointsLayerConfig, PointsLayerUpdate } from "@/react/spatialdata/points_layer_config";

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
    config: PointsLayerConfig;
    updateLayer: PointsLayerUpdate;
};

type CatalogEntry = { code: number; name: string; count?: number };

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

export default function PointsFeatureFilterPanel({ config, updateLayer }: Props) {
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
    // worker, and `loadPointsMatchingFeatureCodes` throws outright without one. MDV
    // cannot enable that worker — the published `@spatialdata/core/points-worker` is a
    // CommonJS file in an ESM package, so `new Worker(url, {type:"module"})` dies on
    // `require is not defined` (SpatialData.js#148). Until that ships fixed, a
    // non-resident feature genuinely cannot be fetched, so say so rather than inviting
    // a click that silently does nothing. Drops out on its own once the worker loads.
    const canScanOnDemand = supportsOnDemandLoad && isPointsWorkerEnabled();
    /** The element could scan, but the worker it needs is unavailable. */
    const workerBlocksScan = supportsOnDemandLoad && !canScanOnDemand;
    const rowInfo = (code: number) => {
        const resident = residentKnown && (residentCodes?.has(code) ?? false);
        const rendered = loadedMatchingCodes?.has(code) ?? false;
        const selected = !noneSelected && (allSelected || selectedCodes.has(code));
        return {
            resident,
            rendered,
            selected,
            state: describeFeatureRowState({
                resident,
                rendered,
                selected,
                scanning,
                supportsOnDemandLoad: canScanOnDemand,
                residentKnown,
            }),
        };
    };
    const notLoadedCount = residentKnown
        ? entries.reduce((total, entry) => total + (rowInfo(entry.code).state.greyed ? 1 : 0), 0)
        : 0;

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
                              ? "not in the loaded sample (greyed below). Fetching them needs the points worker, which this build can't start (SpatialData.js#148) — raise the memory cap to bring more in."
                              : "not in the loaded sample (greyed below). This dataset has no feature index, so they can't be shown until the memory cap is raised or it's rewritten with one."}
                    </Typography>
                </Alert>
            ) : null}

            {matchingLoadState ? (
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
                    const countStr = entry.count !== undefined ? ` · ${entry.count.toLocaleString()} pts` : "";
                    // Multi-line diagnostic: the human state and its reason, then the raw
                    // signals that drove it (what made this row grey, or not).
                    const title =
                        `${entry.name} · code ${entry.code}${countStr}\n` +
                        `${state.label}: ${state.reason}\n` +
                        // The classifier is told there is no on-demand load, which is true
                        // here but for the wrong reason — it blames a missing feature index,
                        // and this element has one. Correct the attribution rather than fork
                        // the classifier, which is a copy of upstream's and due for deletion.
                        (workerBlocksScan
                            ? "(This element does have a feature index; the points worker it needs can't start — SpatialData.js#148.)\n"
                            : "") +
                        `[resident=${resident ? "y" : "n"} rendered=${rendered ? "y" : "n"} ` +
                        `selected=${selected ? "y" : "n"} scan=${scanning ? "running" : "idle"}]`;
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
                                                countIsPartial(entry)
                                                    ? "Points loaded so far (resident window) — dataset total still counting"
                                                    : "Points in the dataset"
                                            }
                                        >
                                            {countIsPartial(entry) ? "≥" : ""}
                                            {effectiveCount(entry)?.toLocaleString() ?? "—"}
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
