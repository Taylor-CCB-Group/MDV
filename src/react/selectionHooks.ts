import { useEffect, useState } from "react";
import { useChart } from "./context";

/**
 * Hook to get the index of the highlighted data point.
 * May want to change how this works in terms of context etc, current implementation is simple
 * will not behave well if used multiple times by the same chart, etc.
 */
export function useHighlightedIndex() {
    const chart = useChart();
    // const { highlightedData } = chart.dataStore;
    const [highlightedIndex, setHighlightedIndex] = useState<number>(-1);

    useEffect(() => {
        console.log("useHighlightedIndex effect");
        // todo check the type - but it's definitely indexable by number, returning number
        chart.onDataHighlighted = ({ indexes }: { indexes: { [k: number]: number } }) => {
            setHighlightedIndex(indexes[0]);
        };
    }, [chart]);

    return highlightedIndex;
}
