import { action } from "mobx";
import { observer } from "mobx-react-lite";
import { useEffect, useRef, useState } from "react";
import { useDataStore } from "../react/context";
import { useChartSize, useConfig } from "../react/hooks";
import { withPhysicalUnit } from "./analysisFrame";
import type { CrossPcfConfig } from "./CrossPcfChart";
import { CrossPcfCurves, type CurvePair, PAIR_COLOURS, PairLegend } from "./CrossPcfCurves";
import { CrossPcfMatrix, matrixAt, type PaintEvent } from "./CrossPcfMatrix";
import { type CrossPcfResult, computeCrossPCF, resolveAnalysisFrame } from "./crossPcf";
import { domainSourceEnv } from "./domainSourceEnv";
import { defaultDomainSourceRef, resolveDomainSourceRef } from "./domainSourceRef";
import { labelGutter, panelLayout } from "./layout";
import { continuePaint, type PaintMode, removePair, resolvePairs, startPaint } from "./pairSelection";
import { addBin, clampBins, moveBin, type RBin, rColour, removeBin, toSelection } from "./rBins";

type State =
    | { status: "idle" | "running" }
    | { status: "done"; result: CrossPcfResult }
    | { status: "error"; message: string };

export const CrossPcfChartComponent = observer(() => {
    const config = useConfig<CrossPcfConfig>();
    const dataStore = useDataStore();
    const [width, height] = useChartSize();

    // Reading the source here (inside `observer`) tracks the other chart's selection and the gate map.
    const env = domainSourceEnv(dataStore);
    const source = resolveDomainSourceRef(config.domain_source ?? defaultDomainSourceRef(env), env);
    const typeColumn = config.type_column || undefined;
    const { r_max: rMax, n_bins: nBins, edge_corrected: edgeCorrected, um_per_unit: umPerUnit } = config;

    const [state, setState] = useState<State>({ status: "idle" });
    // The last finished result stays on screen while a newer one computes (lasso drags edit the source often).
    const [shown, setShown] = useState<CrossPcfResult | null>(null);

    useEffect(() => {
        if (!source) {
            setState({ status: "idle" });
            setShown(null);
            return;
        }
        const abort = new AbortController();
        setState({ status: "running" });
        (async () => {
            const resolved = await resolveAnalysisFrame(dataStore, source.region);
            const frame =
                umPerUnit > 0 && resolved.unit !== "µm"
                    ? withPhysicalUnit(resolved, umPerUnit, "chart setting")
                    : resolved;
            return computeCrossPCF(dataStore, source, {
                typeColumn,
                frame,
                rMax,
                nBins,
                edgeCorrected,
                signal: abort.signal,
            });
        })().then(
            (result) => {
                if (abort.signal.aborted) return;
                setState({ status: "done", result });
                setShown(result);
            },
            (e) => {
                if (!abort.signal.aborted)
                    setState({ status: "error", message: e instanceof Error ? e.message : String(e) });
            },
        );
        return () => abort.abort();
    }, [dataStore, source, typeColumn, rMax, nBins, edgeCorrected, umPerUnit]);

    const res = shown;
    const inside = res ? res.candidates - (res.outside ?? 0) : 0;
    const headerH = 52;

    return (
        <div className="flex flex-col w-full h-full text-sm" style={{ color: "var(--text_color)" }}>
            <div className="px-2 pt-1" style={{ height: headerH }}>
                <div className="truncate">
                    {source ? source.label : "No domain: draw a selection or save a gate, or pick one in settings"}
                    {state.status === "running" && " — computing…"}
                    {state.status === "error" && <span className="text-red-500"> — {state.message}</span>}
                </div>
                {res && (
                    <div className="flex items-center gap-2 opacity-80">
                        <span className="truncate" title={`frame: ${res.frame.provenance}`}>
                            {`${inside} cells, ${res.labels.length} ${res.typeColumn} types, ${res.ms.toFixed(0)} ms`}
                        </span>
                    </div>
                )}
            </div>
            {res && <Panels res={res} width={width} height={Math.max(0, height - headerH)} />}
        </div>
    );
});

/** Curves for the chosen pairs beside the matrix (stacked when the chart is narrow), as in the mock. */
type Pair = { a: number; b: number };

/**
 * Curves (with their legend) beside the matrix — stacked when the chart is narrow — sharing one
 * focused pair: hovering a curve, a legend item or a matrix cell highlights that pair in all three.
 */
const Panels = observer(({ res, width, height }: { res: CrossPcfResult; width: number; height: number }) => {
    const config = useConfig<CrossPcfConfig>();
    const [focus, setFocus] = useState<Pair | null>(null);
    // The r highlighted across the chart: its band, chip and matrix stand out, the others recede.
    const [rFocus, setRFocus] = useState<number | null>(null);
    const stroke = useRef<PaintMode>("add");
    const label = (t: number) => res.labels[t] ?? String(t);
    const asLabels = (p: Pair) => ({ a: label(p.a), b: label(p.b) });
    // Colour follows the pair's position in config, so it doesn't jump when another domain lacks a type.
    const chosen: CurvePair[] = resolvePairs(res.labels, config.pairs).map((p) => ({
        a: p.a,
        b: p.b,
        label: `${label(p.a)} → ${label(p.b)}`,
        colour: PAIR_COLOURS[p.index % PAIR_COLOURS.length],
    }));
    const focusChosen = focus && chosen.some((p) => p.a === focus.a && p.b === focus.b);
    const preview: CurvePair | undefined =
        focus && !focusChosen
            ? { ...focus, label: `${label(focus.a)} → ${label(focus.b)}`, colour: "currentColor" }
            : undefined;
    const paint = action((e: PaintEvent) => {
        if (e.kind === "start") {
            const s = startPaint(config.pairs, asLabels(e.pair), e.modifier);
            stroke.current = s.mode;
            config.pairs = s.pairs;
        } else {
            const next = continuePaint(config.pairs, asLabels(e.pair), stroke.current);
            // moving over a cell the stroke already covers changes nothing: don't churn the config
            if (next.length !== config.pairs.length) config.pairs = next;
        }
    });
    const remove = action((p: Pair) => {
        config.pairs = removePair(config.pairs, asLabels(p));
        setFocus(null);
    });

    // r's: one band on the curves and one matrix each
    const bins = clampBins(config.r_bins, res.nBins);
    const rs = bins.map((bin, i) => ({ bin, colour: rColour(i) }));
    const moveR = action((index: number, k: number) => {
        config.r_bins = moveBin(config.r_bins, index, k);
    });
    const addR = action((b: RBin) => {
        config.r_bins = addBin(config.r_bins, b);
    });
    const removeR = action((index: number) => {
        config.r_bins = removeBin(config.r_bins, index);
        setRFocus(null); // indices shift
    });
    const rText = (b: RBin) =>
        b === "mean"
            ? `mean, r < ${(res.r[res.nBins - 1] + res.r[0]).toFixed(0)} ${res.frame.unit}`
            : `r = ${res.r[b].toFixed(1)} ${res.frame.unit}`;

    const barH = 26; // r chips
    const readoutH = 22;
    const bodyH = Math.max(0, height - barH - readoutH);
    const m = rs.length;
    const view = matrixAt(res, { kind: "mean" }); // only for which types have cells (same for every r)
    const n = view.types.length;
    const labelSize = labelGutter(Math.max(1, ...view.types.map((t) => label(t).length)), 12);
    const L = panelLayout(width, bodyH, m, n, labelSize);
    const { cols, cell, gap, strip } = L.grid;
    const matrixPx = n * cell;

    // One readout for the whole chart: the focused pair's g at every r shown, else a hint.
    const readout = (() => {
        if (!focus)
            return "row = A (anchor), column = B (neighbour) · drag over cells to paint pairs, Shift/⌘ to add or remove";
        const at = (b: RBin) => {
            const v = matrixAt(res, toSelection(b, res.nBins));
            const i = v.types.indexOf(focus.a);
            const j = v.types.indexOf(focus.b);
            if (i < 0 || j < 0) return "—";
            const k = i * v.types.length + j;
            return `g ${v.value[k].toFixed(2)} (${v.pairs[k]} pairs) at ${b === "mean" ? "mean" : `${res.r[b].toFixed(1)} ${res.frame.unit}`}`;
        };
        return `${label(focus.a)} → ${label(focus.b)}: ${bins.map(at).join(" · ")}`;
    })();

    const chipButton = { background: "none", border: "none", color: "inherit", cursor: "pointer" } as const;
    return (
        <div className="flex flex-col text-sm">
            <div className="flex items-center gap-2 px-2" style={{ height: barH }}>
                {rs.map(({ bin, colour }, i) => (
                    <span
                        key={`chip-${bin}`}
                        className="flex items-center gap-1 px-1.5 rounded whitespace-nowrap"
                        style={{
                            border: `${rFocus === i ? 2 : 1}px solid ${colour}`,
                            background: rFocus === i ? `${colour}33` : undefined,
                            fontWeight: rFocus === i ? 700 : 400,
                            opacity: rFocus !== null && rFocus !== i ? 0.6 : 1,
                        }}
                        onPointerEnter={() => setRFocus(i)}
                        onPointerLeave={() => setRFocus(null)}
                    >
                        <span aria-hidden style={{ width: 9, height: 9, background: colour, borderRadius: 2 }} />
                        {rText(bin)}
                        {m > 1 && (
                            <button
                                type="button"
                                aria-label={`Remove ${rText(bin)}`}
                                className="opacity-60 hover:opacity-100"
                                style={chipButton}
                                onClick={() => removeR(i)}
                            >
                                ×
                            </button>
                        )}
                    </span>
                ))}
                {!bins.includes("mean") && (
                    <button
                        type="button"
                        className="opacity-70 hover:opacity-100 px-1.5 rounded whitespace-nowrap"
                        style={{ ...chipButton, border: "1px dashed currentColor" }}
                        onClick={() => addR("mean")}
                    >
                        + mean over r
                    </button>
                )}
                <span className="opacity-60 truncate text-xs">
                    drag a band on the curves to move it, Shift-press to add one
                </span>
            </div>
            <div className="px-2 truncate opacity-90" style={{ height: readoutH }} title={readout}>
                {readout}
            </div>
            {/* scrolls only when the matrices hit their minimum cell size */}
            <div
                className={L.stacked ? "flex flex-col" : "flex flex-row"}
                style={{ gap, height: bodyH, overflow: "auto" }}
            >
                <div className="flex flex-row" style={{ width: L.curves.width, height: L.curves.height }}>
                    <CrossPcfCurves
                        result={res}
                        pairs={chosen}
                        preview={preview}
                        focus={focus}
                        onFocus={setFocus}
                        rs={rs}
                        onMoveR={moveR}
                        onAddR={addR}
                        rFocus={rFocus}
                        onRFocus={setRFocus}
                        unit={res.frame.unit}
                        width={L.curves.width - L.legendWidth}
                        height={L.curves.height}
                    />
                    <PairLegend
                        result={res}
                        pairs={chosen}
                        preview={preview}
                        focus={focus}
                        onFocus={setFocus}
                        onRemove={remove}
                        height={L.curves.height}
                    />
                </div>
                {/* Matrices share one label gutter: row labels left of the first column, column labels
                    above the first row, so every cell lines up across the grid. */}
                <div
                    className="grid"
                    style={{
                        gridTemplateColumns: [labelSize + matrixPx, ...Array(Math.max(0, cols - 1)).fill(matrixPx)]
                            .map((w) => `${w}px`)
                            .join(" "),
                        columnGap: gap,
                        rowGap: gap,
                        alignContent: "start",
                    }}
                >
                    {rs.map(({ bin, colour }, i) => {
                        const firstCol = i % cols === 0;
                        const firstRow = i < cols;
                        return (
                            <div
                                key={`matrix-${bin}`}
                                onPointerEnter={() => setRFocus(i)}
                                onPointerLeave={() => setRFocus(null)}
                                style={{
                                    // the other r's recede a little while one is in focus
                                    opacity: rFocus !== null && rFocus !== i ? 0.55 : 1,
                                    transition: "opacity 120ms",
                                }}
                            >
                                <div
                                    className="flex items-center gap-1 px-1"
                                    style={{
                                        height: strip,
                                        marginLeft: firstCol ? labelSize : 0,
                                        borderTop: `${rFocus === i ? 5 : 3}px solid ${colour}`,
                                        background: rFocus === i ? `${colour}26` : undefined,
                                        fontWeight: rFocus === i ? 700 : 400,
                                    }}
                                >
                                    <span className="truncate flex-1">{rText(bin)}</span>
                                    {m > 1 && (
                                        <button
                                            type="button"
                                            aria-label={`Remove matrix at ${rText(bin)}`}
                                            className="opacity-60 hover:opacity-100"
                                            style={chipButton}
                                            onClick={() => removeR(i)}
                                        >
                                            ×
                                        </button>
                                    )}
                                </div>
                                <CrossPcfMatrix
                                    result={res}
                                    labels={res.labels}
                                    r={toSelection(bin, res.nBins)}
                                    cell={cell}
                                    labelSize={labelSize}
                                    rowLabels={firstCol}
                                    colLabels={firstRow}
                                    selected={chosen}
                                    focus={focus}
                                    onFocus={setFocus}
                                    onPaint={paint}
                                />
                            </div>
                        );
                    })}
                </div>
            </div>
        </div>
    );
});
