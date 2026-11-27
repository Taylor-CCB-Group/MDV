/**
 * PROTOTYPE (layout lab) — throwaway. One view holding every chart, laid out by gridstack 14
 * or by the in-house engine, with a bar to switch engines, pin a width, save and reload the
 * layout in the v2 shape (ADR 0026) and read measurements.
 */
import { type JSX, useCallback, useEffect, useMemo, useState, useSyncExternalStore } from "react";
import type ChartManager from "@/charts/ChartManager";
import type DataStore from "@/datastore/DataStore";
import { type Layout12, sameLayout } from "./engine";
import GridstackMode from "./GridstackMode";
import InHouseMode from "./InHouseMode";
import { clearSaved, collectCharts, convertPaneView, loadSaved, type Measurement, save } from "./labShared";
import "./layoutLab.css";

export type LabMode = "inhouse" | "gridstack";
const ENGINE_NAMES: Record<LabMode, string> = {
    inhouse: "In-house",
    gridstack: "gridstack 14",
};
const WIDTHS = ["fit", "1440", "1024", "768", "390"] as const;

let nextId = 0;
function useFilterSize(dataStore: DataStore) {
    const subscribe = useCallback(
        (onChange: () => void) => {
            const id = `layout-lab-${nextId++}`;
            dataStore.addListener(id, (type: string) => {
                if (type === "filtered") onChange();
            });
            return () => dataStore.removeListener(id);
        },
        [dataStore],
    );
    return useSyncExternalStore(subscribe, () => dataStore.filterSize);
}

function DataSourceChip({ name, dataStore }: { name: string; dataStore: DataStore }) {
    const filtered = useFilterSize(dataStore);
    return (
        <span className="lab-chip" title={`${name}: filtered rows / total`}>
            {name} <b>{filtered.toLocaleString()}</b> / {dataStore.size.toLocaleString()}
        </span>
    );
}

/** Bumps whenever the view or its charts change, so the lab re-collects them. */
function useViewVersion(cm: ChartManager) {
    const [version, setVersion] = useState(0);
    useEffect(() => {
        const id = `layout-lab-${nextId++}`;
        cm.addListener(id, (type: string) => {
            if (type === "view_loaded" || type === "chart_added" || type === "chart_removed") setVersion((v) => v + 1);
        });
        return () => cm.removeListener(id);
    }, [cm]);
    return version;
}

type LabWindow = {
    __labLayout?: Layout12;
    __labSaved?: Layout12 | null;
    __labConverted?: Layout12;
};

export default function LayoutLab({ cm, initialMode }: { cm: ChartManager; initialMode: LabMode }) {
    const [mode, setMode] = useState<LabMode>(initialMode);
    const [width, setWidth] = useState<(typeof WIDTHS)[number]>("fit");
    const [results, setResults] = useState<Measurement[]>([]);
    const version = useViewVersion(cm);
    // remounts the engine when the layout is replaced from outside it
    const [gen, setGen] = useState(0);

    // biome-ignore lint/correctness/useExhaustiveDependencies: re-collect on view changes
    const charts = useMemo(() => collectCharts(cm), [cm, version]);
    const converted = useMemo(() => convertPaneView(cm, charts), [cm, charts]);
    // biome-ignore lint/correctness/useExhaustiveDependencies: a new view starts from what it saved
    const [saved, setSaved] = useState<Layout12 | null>(() => loadSaved(cm));
    const [layout, setLayout] = useState<Layout12>(() => saved ?? converted);
    // the view whose layout is loaded; engines wait for it so they never mount with the last view's layout
    const [loadedVersion, setLoadedVersion] = useState(version);
    // biome-ignore lint/correctness/useExhaustiveDependencies: reload when the view changes
    useEffect(() => {
        const s = loadSaved(cm);
        setSaved(s);
        setLayout(s ?? converted);
        setLoadedVersion(version);
    }, [version]);

    const w = window as unknown as LabWindow;
    w.__labLayout = layout;
    w.__labSaved = saved;
    w.__labConverted = converted;

    const dirty = !saved || !sameLayout(saved, layout);
    const onMeasured = (m: Measurement) => setResults((r) => [m, ...r].slice(0, 4));
    const props = { charts, layout, onLayoutChange: setLayout, onMeasured };
    const layoutKey = `${mode}-${version}-${gen}`;
    const engine = (
        {
            inhouse: <InHouseMode key={layoutKey} {...props} />,
            gridstack: <GridstackMode key={layoutKey} {...props} />,
        } satisfies Record<LabMode, JSX.Element>
    )[mode];
    const last = results[0];
    return (
        <div className="lab-root">
            <div className="lab-strip" role="toolbar" aria-label="Layout lab">
                <span className="lab-strip-title">One view · {charts.length} charts</span>
                {cm.dataSources.map((ds) => (
                    <DataSourceChip key={ds.name} name={ds.name} dataStore={ds.dataStore} />
                ))}
                <span className="lab-spacer" />
                <span className="lab-tag">Layout lab</span>
                <div className="lab-seg">
                    {(Object.keys(ENGINE_NAMES) as LabMode[]).map((m) => (
                        <button key={m} type="button" data-mode={m} aria-pressed={mode === m} onClick={() => setMode(m)}>
                            {ENGINE_NAMES[m]}
                        </button>
                    ))}
                </div>
                <div className="lab-seg">
                    {WIDTHS.map((x) => (
                        <button key={x} type="button" aria-pressed={width === x} onClick={() => setWidth(x)}>
                            {x === "fit" ? "Fit" : `${x}px`}
                        </button>
                    ))}
                </div>
                <button
                    type="button"
                    className="lab-act"
                    data-act="save"
                    onClick={() => {
                        save(cm, layout);
                        setSaved(layout);
                    }}
                >
                    Save{dirty && <span className="lab-dirty"> •</span>}
                </button>
                <button
                    type="button"
                    className="lab-act"
                    data-act="reset"
                    onClick={() => {
                        clearSaved(cm);
                        setSaved(null);
                        setLayout(converted);
                        setGen((g) => g + 1);
                    }}
                >
                    Reconvert
                </button>
                <button
                    type="button"
                    className="lab-act"
                    data-act="copy"
                    onClick={() => void navigator.clipboard?.writeText(JSON.stringify({ version: 2, layout }, null, 2))}
                >
                    Copy layout JSON
                </button>
                <span className="lab-results" title="Last drag or resize">
                    {last
                        ? `${last.avgFps} fps · worst ${last.worstFrameMs} ms · ${last.setSizeCalls} redraws`
                        : "Drag or resize a chart to measure"}
                </span>
            </div>
            <div className="lab-stage">
                <div className="lab-frame" style={{ maxWidth: width === "fit" ? "100%" : `${width}px` }}>
                    {loadedVersion === version ? engine : null}
                </div>
            </div>
        </div>
    );
}
