/**
 * PROTOTYPE (layout lab) — throwaway. gridstack 14 (aliased as `gridstack14`, so MDV's own
 * gridstack 7.3 is untouched) hosting real charts. 12 columns from 960 px, 6 from 600 px and 1
 * below. Narrower layouts come from the prototype's own `flow`, passed as gridstack's column
 * layout callback (its built-in `list` mode made every chart full width at 6 columns). Moving
 * and resizing only at 12 columns.
 */
import { GridStack, type GridStackWidget } from "gridstack14";
import "gridstack14/dist/gridstack.min.css";
import { useEffect, useLayoutEffect, useRef } from "react";
import ChartHost from "./ChartHost";
import { COLS, flow, GAP_PX, type Layout12, ROW_PX } from "./engine";
import { gestureStart, gestureStop, type LabChart, type Measurement, meterStart, meterStop } from "./labShared";

type Props = {
    charts: LabChart[];
    layout: Layout12;
    onLayoutChange: (layout: Layout12) => void;
    onMeasured: (m: Measurement) => void;
};

export default function GridstackMode({ charts, layout, onLayoutChange, onMeasured }: Props) {
    const gridRef = useRef<HTMLDivElement>(null);
    const grid = useRef<GridStack | null>(null);
    const callbacks = useRef({ onLayoutChange, onMeasured });
    callbacks.current = { onLayoutChange, onMeasured };
    // the layout is read once, when gridstack takes over the items; afterwards gridstack owns it
    const initial = useRef(layout);
    // the last 12-column layout, which narrower layouts are flowed from
    const wide = useRef(layout);

    useLayoutEffect(() => {
        const el = gridRef.current;
        if (!el) return;
        const g = GridStack.init(
            {
                column: COLS,
                cellHeight: ROW_PX + GAP_PX,
                margin: GAP_PX / 2,
                animate: true,
                handle: ".ciview-chart-titlebar",
                draggable: { cancel: ".ciview-chart-icon,.ciview-menu-icon,input,textarea,button,select,option" },
                resizable: { handles: "e,se,s,sw,w", autoHide: false },
                columnOpts: {
                    columnMax: COLS,
                    // gridstack hands over the items still to place (`nodes`) and expects them pushed
                    // into `newNodes`, which already holds any it restored from its own cache
                    layout: (column, _old, newNodes, nodes) => {
                        const placed = column === COLS ? wide.current : null;
                        const flowed = column === COLS ? null : flow(wide.current, column);
                        for (const n of nodes) {
                            const id = String(n.id);
                            const r = placed ? placed[id] : flowed?.find((o) => o.i === id);
                            if (r) {
                                n.x = r.x;
                                n.y = r.y;
                                n.w = r.w;
                                n.h = r.h;
                            }
                            newNodes.push(n);
                        }
                    },
                    breakpoints: [
                        { w: 599, c: 1 },
                        { w: 959, c: 6 },
                    ],
                },
            },
            el,
        );
        if (!g) return;
        grid.current = g;
        const editable = () => g.getColumn() === COLS;
        const applyEditable = () => {
            const on = editable();
            g.enableMove(on);
            g.enableResize(on);
            (window as unknown as { __labCols?: number }).__labCols = g.getColumn();
        };
        applyEditable();
        g.on("change", () => {
            if (!editable()) {
                applyEditable();
                return;
            }
            applyEditable();
            const saved = g.save(false) as GridStackWidget[];
            const next: Layout12 = {};
            for (const w of saved) {
                if (w.id == null) continue;
                next[String(w.id)] = { x: w.x ?? 0, y: w.y ?? 0, w: w.w ?? 1, h: w.h ?? 1 };
            }
            wide.current = next;
            callbacks.current.onLayoutChange(next);
        });
        const start = () => {
            gestureStart();
            meterStart();
        };
        const stop = (kind: string) => () => {
            gestureStop();
            void meterStop(`gridstack ${kind} @ ${g.getColumn()} cols`).then(callbacks.current.onMeasured);
        };
        g.on("dragstart", start);
        g.on("resizestart", start);
        g.on("dragstop", stop("move"));
        g.on("resizestop", stop("resize"));
        return () => {
            g.destroy(false);
            grid.current = null;
        };
    }, []);

    // a column change from the width presets fires no `change` when nothing moved
    useEffect(() => {
        const el = gridRef.current;
        if (!el) return;
        const ro = new ResizeObserver(() => {
            const g = grid.current;
            if (!g) return;
            const on = g.getColumn() === COLS;
            g.enableMove(on);
            g.enableResize(on);
            (window as unknown as { __labCols?: number }).__labCols = g.getColumn();
        });
        ro.observe(el);
        return () => ro.disconnect();
    }, []);

    return (
        <div className="lab-scroll lab-canvas">
            <div ref={gridRef} className="grid-stack lab-gs">
                {charts.map((c) => {
                    const r = initial.current[c.id] ?? { x: 0, y: 1000, w: 6, h: 6 };
                    return (
                        <div
                            key={c.id}
                            className="grid-stack-item"
                            data-id={c.id}
                            gs-id={c.id}
                            gs-x={r.x}
                            gs-y={r.y}
                            gs-w={r.w}
                            gs-h={r.h}
                        >
                            <div className="grid-stack-item-content lab-card">
                                <ChartHost chart={c.chart} />
                            </div>
                        </div>
                    );
                })}
            </div>
        </div>
    );
}
