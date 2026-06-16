import { observer } from "mobx-react-lite";
import { Fragment, useEffect, useState } from "react";
import { v4 as uuidv4 } from "uuid";
import { useChart } from "../../react/context";
import { useChartManager } from "../../react/hooks";

const MAX_VISIBLE_CELLS = 2000;

function useFilter(datastore) {
    const [filter, setFilter] = useState({});

    useEffect(() => {
        const id = uuidv4();
        datastore.addListener(id, (type, data) => {
            if (type === "filtered") {
                setFilter({
                    timeStamp: Date.now(),
                    dimension: data,
                });
            }
        });

        return () => {
            datastore.removeListener(id);
        };
    }, [datastore]);

    return filter;
}


// Blue color scale for nonzero, neutral gray for zero
function getCellColor(value, maxValue, isDark) {
    if (maxValue <= 0 || value <= 0) {
        return isDark ? "#23272e" : "#e0e2e6";
    }
    // Use a blue scale (light to dark)
    const intensity = value / maxValue;
    if (isDark) {
        // Blue: low = #2a3556, high = #4fc3f7
        const low = [42, 53, 86]; // #2a3556
        const high = [79, 195, 247]; // #4fc3f7
        const rgb = low.map((l, i) => Math.round(l + (high[i] - l) * intensity));
        return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
    } else {
        // Blue: low = #e0e8f7, high = #1976d2
        const low = [224, 232, 247]; // #e0e8f7
        const high = [25, 118, 210]; // #1976d2
        const rgb = low.map((l, i) => Math.round(l + (high[i] - l) * intensity));
        return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
    }
}

const cellBaseStyle = {
    minWidth: "48px",
    minHeight: "32px",
    border: "1px solid var(--heatmap-border, #d5dde5)",
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    fontSize: "12px",
    transition: "background 0.2s, color 0.2s",
};

function pruneAggregation(aggregation, xDisplayCategories, yDisplayCategories) {
    const normalizeSelection = (selected) =>
        (selected || [])
            .map((value) => {
                if (typeof value === "string") {
                    return value;
                }
                if (value && typeof value === "object") {
                    if (typeof value.value === "string") {
                        return value.value;
                    }
                    if (typeof value.label === "string") {
                        return value.label;
                    }
                }
                return null;
            })
            .filter((value) => typeof value === "string");

    const xFilter = new Set(normalizeSelection(xDisplayCategories));
    const yFilter = new Set(normalizeSelection(yDisplayCategories));
    const useXFilter = xFilter.size > 0;
    const useYFilter = yFilter.size > 0;

    const xIndexes = [];
    const xLabels = [];
    for (let i = 0; i < aggregation.xLabels.length; i++) {
        const label = aggregation.xLabels[i];
        if (!useXFilter || xFilter.has(label)) {
            xIndexes.push(i);
            xLabels.push(label);
        }
    }

    const yIndexes = [];
    const yLabels = [];
    for (let i = 0; i < aggregation.yLabels.length; i++) {
        const label = aggregation.yLabels[i];
        if (!useYFilter || yFilter.has(label)) {
            yIndexes.push(i);
            yLabels.push(label);
        }
    }

    const counts = [];
    let maxCount = 0;
    for (const xIndex of xIndexes) {
        const srcRow = aggregation.counts[xIndex] || [];
        const row = [];
        for (const yIndex of yIndexes) {
            const value = srcRow[yIndex] || 0;
            row.push(value);
            if (value > maxCount) {
                maxCount = value;
            }
        }
        counts.push(row);
    }

    return {
        xLabels,
        yLabels,
        counts,
        maxCount,
        totalCells: xLabels.length * yLabels.length,
    };
}

const CategoryHeatmapComponent = observer(() => {
    const chartManager = useChartManager();
    const isDark = chartManager?.theme === "dark";
    const chart = useChart();
    const filter = useFilter(chart.dataStore);
    const paramKey = JSON.stringify(chart.config.param || []);
    const [aggregation, setAggregation] = useState(null);
    const [isLoading, setIsLoading] = useState(true);
    const [errorMessage, setErrorMessage] = useState("");

    useEffect(() => {
        let cancelled = false;
        setIsLoading(true);
        setErrorMessage("");

        chart
            .fetchAggregation()
            .then((nextAggregation) => {
                if (cancelled) {
                    return;
                }
                setAggregation(nextAggregation);
            })
            .catch((error) => {
                if (cancelled) {
                    return;
                }
                setAggregation(null);
                setErrorMessage(error?.message || "Failed to aggregate category heatmap data.");
            })
            .finally(() => {
                if (!cancelled) {
                    setIsLoading(false);
                }
            });

        return () => {
            cancelled = true;
        };
    }, [chart, filter.timeStamp, paramKey]);

    const displayData =
        aggregation == null
            ? null
            : pruneAggregation(
                  aggregation,
                  chart.config.x_display_categories,
                  chart.config.y_display_categories,
              );

    if (isLoading) {
        return <div style={{ padding: "12px" }}>Loading category heatmap...</div>;
    }

    if (errorMessage) {
        return (
            <div style={{ padding: "12px", color: "#b00020" }}>{errorMessage}</div>
        );
    }

    if (!displayData || displayData.xLabels.length === 0 || displayData.yLabels.length === 0) {
        return <div style={{ padding: "12px" }}>No category pairs available to display.</div>;
    }

    if (displayData.totalCells > MAX_VISIBLE_CELLS) {
        return (
            <div style={{ padding: "12px", color: "#8a1c0f" }}>
                Category heatmap not drawn: {displayData.yLabels.length} x {displayData.xLabels.length} = {displayData.totalCells} cells exceeds the limit of {MAX_VISIBLE_CELLS}.
            </div>
        );
    }

    // Axis swap: yLabels as rows, xLabels as columns
    return (
        <div style={{ width: "100%", height: "100%" }}>
            <div
                style={{
                    display: "grid",
                    gridTemplateColumns: `minmax(140px, auto) repeat(${displayData.xLabels.length}, minmax(48px, 1fr))`,
                    gap: "0px",
                    overflow: "auto",
                    width: "100%",
                    height: "calc(100% - 40px)",
                    padding: "8px",
                    boxSizing: "border-box",
                    background: isDark ? "#181a20" : undefined,
                }}
            >
                {/* Top-left corner cell */}
                <div
                    style={{
                        ...cellBaseStyle,
                        fontWeight: 600,
                        background: isDark ? "#23272e" : "#f0f3f7",
                        color: isDark ? "#e0e6ef" : undefined,
                    }}
                >
                    Category
                </div>
                {/* Column headers (xLabels) */}
                {displayData.xLabels.map((label) => (
                    <div
                        key={`x-header-${label}`}
                        style={{
                            ...cellBaseStyle,
                            fontWeight: 600,
                            background: isDark ? "#23272e" : "#f0f3f7",
                            color: isDark ? "#e0e6ef" : undefined,
                            writingMode: "vertical-rl",
                            transform: "rotate(180deg)",
                        }}
                        title={label}
                    >
                        {label}
                    </div>
                ))}
                {/* Data rows: yLabels as rows */}
                {displayData.yLabels.map((yLabel, yIndex) => (
                    <Fragment key={`y-row-${yLabel}`}>
                        {/* Row header */}
                        <div
                            key={`y-label-${yLabel}`}
                            style={{
                                ...cellBaseStyle,
                                justifyContent: "flex-start",
                                padding: "0 8px",
                                background: isDark ? "#23272e" : "#fafbfd",
                                color: isDark ? "#e0e6ef" : undefined,
                                fontWeight: 500,
                            }}
                            title={yLabel}
                        >
                            {yLabel}
                        </div>
                        {/* Data cells */}
                        {displayData.xLabels.map((xLabel, xIndex) => {
                            const value = displayData.counts[xIndex]?.[yIndex] || 0;
                            const isActive =
                                chart.cellFilter &&
                                chart.cellFilter[0] === xLabel &&
                                chart.cellFilter[1] === yLabel;
                            return (
                                <div
                                    key={`cell-${xLabel}-${yLabel}`}
                                    style={{
                                        ...cellBaseStyle,
                                        background: getCellColor(value, displayData.maxCount, isDark),
                                        color: value > 0 ? (isDark ? "#e0e6ef" : "#0e2238") : (isDark ? "#7a8699" : "#8a96a3"),
                                        cursor: "pointer",
                                        boxShadow: isActive
                                            ? (isDark ? "inset 0 0 0 2px #4fc3f7" : "inset 0 0 0 2px #0b4f8a")
                                            : "none",
                                    }}
                                    title={`${xLabel} x ${yLabel}: ${value}`}
                                    onClick={() => chart.filterCell(xLabel, yLabel)}
                                >
                                    {value}
                                </div>
                            );
                        })}
                    </Fragment>
                ))}
            </div>
            {/* Legend */}
            <div style={{ display: "flex", alignItems: "center", gap: 8, marginTop: 8, paddingLeft: 8 }}>
                <span style={{ fontSize: 12, color: isDark ? "#e0e6ef" : "#0e2238" }}>Low</span>
                <div style={{ display: "flex", height: 16, width: 120, background: "none" }}>
                    {Array.from({ length: 20 }).map((_, i) => {
                        const v = i / 19;
                        return (
                            <div
                                key={i}
                                style={{
                                    flex: 1,
                                    background: getCellColor(1 + v * 99, 100, isDark),
                                    height: "100%",
                                }}
                            />
                        );
                    })}
                </div>
                <span style={{ fontSize: 12, color: isDark ? "#e0e6ef" : "#0e2238" }}>High</span>
                <span style={{ fontSize: 12, color: isDark ? "#8a96a3" : "#7a8699" }}>&nbsp;&nbsp;0</span>
                <div style={{ width: 16, height: 16, background: getCellColor(0, 100, isDark), border: "1px solid #aaa", display: "inline-block", marginLeft: 2 }} />
            </div>
        </div>
    );
});

export default CategoryHeatmapComponent;
