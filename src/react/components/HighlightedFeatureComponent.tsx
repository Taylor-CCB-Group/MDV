import { useEffect, useState } from "react";
import { useChart, useDataStore } from "../context";
import { observer } from "mobx-react-lite";
import type { ColumnName } from "../../charts/charts";
import { useHighlightedIndex } from "../selectionHooks";

function useColumnData(columnName: ColumnName, maxItems = 100) {
    const ds = useDataStore();
    const [columnData, setColumnData] = useState<Float32Array>(null);
    const [indices, setIndices] = useState<Uint32Array>(null);
    const [column, setColumn] = useState(null);
    // biome-ignore lint/correctness/useExhaustiveDependencies: ds._filteredIndicesPromise is used to trigger reactivity
    useEffect(() => {
        if (!window.mdv.chartManager) return;
        const cm = window.mdv.chartManager;
        const indexPromise = ds.getFilteredIndices();
        const colPromise = cm._getColumnsAsync(ds.name, [columnName]);
        setColumn(null);
        setColumnData(null);
        setIndices(null);
        let cancel = false;
        console.log("useColumnData effect");
        Promise.all([indexPromise, colPromise]).then(([indices]) => {
            console.log("useColumnData promise resolved");
            if (cancel) return;
            const column = ds.columnIndex[columnName];
            setColumn(column);
            const newColumnData = new Float32Array(indices.length);
            for (let j = 0; j < indices.length; j++) {
                const i = indices[j];
                const v = (newColumnData[j] = column.data[i]);
                if (Number.isNaN(v)) {
                    newColumnData[j] = Number.NEGATIVE_INFINITY;
                }
            }
            // may want to have a JS array of (index, value) tuples instead of two typed arrays...
            // we end up spreading the columnData at the moment, potentially multiple times...
            // and duplicate work in sorting the indices...
            const sortedIndices = indices.sort(
                (a, b) => column.data[b] - column.data[a],
            );
            setIndices(sortedIndices.slice(0, maxItems));
            setColumnData(
                newColumnData.sort((a, b) => b - a).slice(0, maxItems),
            );
        });
        return () => {
            cancel = true;
        };
    }, [
        columnName,
        ds._filteredIndicesPromise,
        ds.name,
        ds.columnIndex,
        maxItems,
        ds.getFilteredIndices,
    ]);
    return { columnData, indices, column };
}

function useMarkdownText(md: string) {
    const [renderTextFn, setRenderTextFn] =
        useState<(md: string) => string>(null);
    const [html, setHtml] = useState<string>("");
    useEffect(() => {
        import("../../utilities/MarkdownText").then(
            ({ default: renderTextFn }) => {
                setRenderTextFn(() => renderTextFn);
            },
        );
    }, []);
    useEffect(() => {
        if (!md) return setHtml("");
        if (renderTextFn) {
            const html = renderTextFn(md);
            setHtml(html);
        }
    }, [md, renderTextFn]);
    return html;
}

export const HighlightedFeatureComponent = observer(() => {
    const chart = useChart();
    const html = useMarkdownText(chart.config.text);
    const { columnData, indices, column } = useColumnData(
        chart.config.param[0],
        8,
    );
    const highlightedIndex = useHighlightedIndex();
    const highlightInRange =
        indices && highlightedIndex >= 0 && indices.includes(highlightedIndex);
    const highlightValue = column?.data[highlightedIndex];
    return (
        <div style={{ fontSize: "0.8rem" }}>
            {/* biome-ignore lint/security/noDangerouslySetInnerHtml: todo use ReactMarkdown component instead */}
            <div dangerouslySetInnerHTML={{ __html: html }} />
            <em>"{chart.config.param[0]}"</em>
            <br />
            {highlightInRange
                ? `Highlighted value (one of the top values): ${highlightValue}`
                : `Highlighted value ${highlightValue} not in top ${columnData?.length} values.`}
            {columnData &&
                [...columnData].map((d, i) => (
                    <div
                        style={{
                            cursor: "pointer",
                            border: `1px ${highlightedIndex === indices[i] ? "solid" : "dashed"} #ccc`,
                            padding: "0.5em",
                            margin: "0.5em",
                            borderRadius: "0.5em",
                        }}
                        onClick={() => {
                            chart.dataStore.dataHighlighted(
                                [indices[i]],
                                chart,
                            );
                        }}
                        key={indices[i]}
                    >
                        {d}
                    </div>
                ))}
        </div>
    );
});
