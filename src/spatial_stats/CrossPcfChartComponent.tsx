import { action } from "mobx";
import { observer } from "mobx-react-lite";
import { useEffect, useRef, useState } from "react";
import { useDataStore } from "../react/context";
import { useChartSize, useConfig } from "../react/hooks";
import { withPhysicalUnit } from "./analysisFrame";
import type { CrossPcfConfig } from "./CrossPcfChart";
import { CrossPcfCurves, type CurvePair, PAIR_COLOURS, PairLegend } from "./CrossPcfCurves";
import { CrossPcfMatrix, type PaintEvent, type RSelection } from "./CrossPcfMatrix";
import { type CrossPcfResult, computeCrossPCF, resolveAnalysisFrame } from "./crossPcf";
import { domainSourceEnv } from "./domainSourceEnv";
import { defaultDomainSourceRef, resolveDomainSourceRef } from "./domainSourceRef";
import { continuePaint, type PaintMode, removePair, resolvePairs, startPaint } from "./pairSelection";

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

    const r: RSelection = config.r_bin === "mean" ? { kind: "mean" } : { kind: "bin", k: config.r_bin };
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
                        <input
                            type="range"
                            aria-label="r bin"
                            min={-1}
                            max={res.nBins - 1}
                            value={config.r_bin === "mean" ? -1 : config.r_bin}
                            onChange={action((e) => {
                                const k = Number(e.currentTarget.value);
                                config.r_bin = k < 0 ? "mean" : k;
                            })}
                        />
                        <span className="whitespace-nowrap">
                            {r.kind === "mean"
                                ? `mean, r < ${(res.r[res.nBins - 1] + res.r[0]).toFixed(0)} ${res.frame.unit}`
                                : `r = ${res.r[r.kind === "bin" ? Math.min(r.k, res.nBins - 1) : 0].toFixed(1)} ${res.frame.unit}`}
                        </span>
                    </div>
                )}
            </div>
            {res && <Panels res={res} r={clampR(r, res.nBins)} width={width} height={Math.max(0, height - headerH)} />}
        </div>
    );
});

const clampR = (r: RSelection, nBins: number): RSelection =>
    r.kind === "bin" ? { kind: "bin", k: Math.min(r.k, nBins - 1) } : r;

/** Curves for the chosen pairs beside the matrix (stacked when the chart is narrow), as in the mock. */
type Pair = { a: number; b: number };

/**
 * Curves (with their legend) beside the matrix — stacked when the chart is narrow — sharing one
 * focused pair: hovering a curve, a legend item or a matrix cell highlights that pair in all three.
 */
const Panels = observer(
    ({ res, r, width, height }: { res: CrossPcfResult; r: RSelection; width: number; height: number }) => {
        const config = useConfig<CrossPcfConfig>();
        const [focus, setFocus] = useState<Pair | null>(null);
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
        const sideBySide = width >= 640;
        const curvesW = sideBySide ? Math.round(width * 0.58) : width;
        const curvesH = sideBySide ? Math.min(height, Math.round(curvesW * 0.62)) : Math.round(height * 0.4);
        const legendW = 120;
        const matrixW = sideBySide ? width - curvesW : width;
        const matrixH = sideBySide ? height : height - curvesH;
        return (
            <div className={sideBySide ? "flex flex-row" : "flex flex-col"}>
                <div className="flex flex-row" style={{ width: curvesW, height: curvesH }}>
                    <CrossPcfCurves
                        result={res}
                        pairs={chosen}
                        preview={preview}
                        focus={focus}
                        onFocus={setFocus}
                        r={r}
                        unit={res.frame.unit}
                        width={curvesW - legendW}
                        height={curvesH}
                    />
                    <PairLegend
                        pairs={chosen}
                        preview={preview}
                        focus={focus}
                        onFocus={setFocus}
                        onRemove={remove}
                        height={curvesH}
                    />
                </div>
                <CrossPcfMatrix
                    result={res}
                    labels={res.labels}
                    r={r}
                    unit={res.frame.unit}
                    width={matrixW}
                    height={matrixH}
                    selected={chosen}
                    focus={focus}
                    onFocus={setFocus}
                    onPaint={paint}
                />
            </div>
        );
    },
);
