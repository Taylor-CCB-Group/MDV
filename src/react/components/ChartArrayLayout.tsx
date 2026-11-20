/**
 * Multi-chart layout shell (density grid and similar). Layout modes, flex/grid CSS, and
 * resize measurement are an initial prototype — not a stable API; expect rework for
 * user-configurable layouts and context-specific arrangements.
 */
import { type CSSProperties, type ReactNode, type Ref, useCallback, useImperativeHandle, useRef } from "react";
import "./chartArrayLayout.css";
import { useChartArrayLayoutGeometry } from "../hooks/useChartArrayLayoutGeometry";
import { type ChartArrayLayoutMode, getChartArrayLayoutMode } from "./chartArrayLayoutUtils";

export type ChartArrayLayoutHandle = {
    root: HTMLDivElement | null;
    cells: readonly (HTMLDivElement | null)[];
};

type ChartArrayLayoutProps = {
    cellCount: number;
    cellKeys?: readonly string[];
    /** Defaults from {@link getChartArrayLayoutMode}; override for future user layout prefs. */
    layout?: ChartArrayLayoutMode;
    layoutRef?: Ref<ChartArrayLayoutHandle>;
    className?: string;
    style?: CSSProperties;
    /** A cell's content, which fills it; a cell's own canvas belongs here. */
    renderCell: (index: number) => ReactNode;
};

export default function ChartArrayLayout({
    cellCount,
    cellKeys,
    layout,
    layoutRef,
    className,
    style,
    renderCell,
}: ChartArrayLayoutProps) {
    const layoutMode = layout ?? getChartArrayLayoutMode(cellCount);
    const rootRef = useRef<HTMLDivElement>(null);
    const { gridFillsViewport, gridColumns } = useChartArrayLayoutGeometry(rootRef, cellCount, layoutMode);
    const cellRefs = useRef<(HTMLDivElement | null)[]>([]);

    const setCellRef = useCallback((index: number, element: HTMLDivElement | null) => {
        cellRefs.current[index] = element;
    }, []);

    useImperativeHandle(
        layoutRef,
        () => ({
            get root() {
                return rootRef.current;
            },
            get cells() {
                return cellRefs.current;
            },
        }),
        [],
    );

    return (
        <div
            ref={rootRef}
            className={["mdv-chart-array", className].filter(Boolean).join(" ")}
            data-layout={layoutMode}
            data-grid-fill={layoutMode === "grid" && gridFillsViewport ? "true" : undefined}
            style={
                layoutMode === "grid"
                    ? {
                          ...style,
                          ["--mdv-chart-array-columns" as string]: String(gridColumns),
                      }
                    : style
            }
        >
            {Array.from({ length: cellCount }, (_, index) => (
                <div
                    key={cellKeys?.[index] ?? index}
                    ref={(element) => setCellRef(index, element)}
                    className="mdv-chart-array__cell"
                    data-chart-array-index={String(index)}
                >
                    {renderCell(index)}
                </div>
            ))}
        </div>
    );
}
