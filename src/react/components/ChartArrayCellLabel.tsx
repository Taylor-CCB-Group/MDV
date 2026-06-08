type ChartArrayCellLabelProps = {
    label: string;
    accentRgb: readonly [number, number, number];
};

export default function ChartArrayCellLabel({ label, accentRgb }: ChartArrayCellLabelProps) {
    return (
        <div
            className="pointer-events-none absolute left-2 top-2 z-10 max-w-[calc(100%-16px)] truncate rounded-sm px-2 py-1 text-xs"
            style={{
                backgroundColor: "rgba(15, 23, 42, 0.76)",
                color: `rgb(${accentRgb[0]}, ${accentRgb[1]}, ${accentRgb[2]})`,
            }}
            title={label}
        >
            {label}
        </div>
    );
}
