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
    renderCell: (index: number) => ReactNode;
};

export default function ChartArrayLayout({
    cellCount,
    cellKeys,
    layoutRef,
    className,
    style,
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
