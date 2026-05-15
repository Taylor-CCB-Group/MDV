import {
    useCallback,
    useImperativeHandle,
    useRef,
    type CSSProperties,
    type ReactNode,
    type Ref,
} from "react";

export type ChartArrayLayoutHandle = {
    root: HTMLDivElement | null;
    cells: readonly (HTMLDivElement | null)[];
};

type ChartArrayLayoutProps = {
    cellCount: number;
    cellKeys?: readonly string[];
    layoutRef?: Ref<ChartArrayLayoutHandle>;
    className?: string;
    style?: CSSProperties;
    /** Shared-canvas layer (e.g. DeckGL) sized to the grid root via CSS. */
    canvasOverlay?: ReactNode;
    renderCell: (index: number) => ReactNode;
};

export default function ChartArrayLayout({
    cellCount,
    cellKeys,
    layoutRef,
    className,
    style,
    canvasOverlay,
    renderCell,
}: ChartArrayLayoutProps) {
    const rootRef = useRef<HTMLDivElement>(null);
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
            style={style}
        >
            {canvasOverlay ? <div className="mdv-chart-array__canvas">{canvasOverlay}</div> : null}
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
