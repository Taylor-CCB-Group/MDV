import { runInAction } from "mobx";
import { useEffect, useId, useMemo } from "react";
import { useDebounce } from "use-debounce";
import type { DataType, LoadedDataColumn } from "../../charts/charts";
import { useChart, useDataStore } from "../context";
import { useDimensionFilter } from "../hooks";
import { featureSearchMatcher, findMissingTerms } from "./featureSearch";
import type { GeneNetworkConfig } from "./GeneNetworkChart";

/** Key the search filter is reported under in `chart.storedFilters`. */
export const FEATURE_SEARCH_FILTER_KEY = "search";

/**
 * Applies the chart's `search` config as a filter on the whole data source, and clears it
 * on the global "remove all filters".
 */
export function useFeatureSearchFilter(column: LoadedDataColumn<DataType>) {
    const chart = useChart<GeneNetworkConfig>();
    const ds = useDataStore();
    const search = chart.config.search;
    // typing narrows the whole view, so let a burst of keystrokes land as one filter
    const [query] = useDebounce(search.query ?? "", 150);
    const terms = search.terms ?? [];
    const termsKey = terms.join("\n");
    const dim = useDimensionFilter(column);

    useEffect(() => {
        if (!dim) return;
        const match = featureSearchMatcher({ query, terms: termsKey ? termsKey.split("\n") : [] });
        if (!match) {
            dim.removeFilter();
            chart.storedFilters.applied(FEATURE_SEARCH_FILTER_KEY);
            return;
        }
        // reads values on each call: the data store reapplies this predicate when the column's data changes
        const predicate = (i: number) => match(String(column.getValue(i)));
        chart.storedFilters.applied(
            FEATURE_SEARCH_FILTER_KEY,
            dim.filter("filterPredicate", [column.field], { predicate }),
        );
    }, [dim, column, query, termsKey, chart]);

    const listenerId = useId();
    useEffect(() => {
        const key = `FeatureSearch-${listenerId}`;
        ds.addListener(key, (type: string, data: unknown) => {
            if (type === "filtered" && data === "all_removed") chart.removeFilter();
        });
        return () => ds.removeListener(key);
    }, [ds, listenerId, chart]);
}

/** Entered terms that are not values of `column`. */
export function useMissingTerms(column: LoadedDataColumn<DataType>, terms: readonly string[]) {
    const ds = useDataStore();
    const termsKey = terms.join("\n");
    return useMemo(() => {
        if (!termsKey) return [];
        function* values() {
            for (let i = 0; i < ds.size; i++) yield String(column.getValue(i));
        }
        return findMissingTerms(termsKey.split("\n"), values());
    }, [ds, column, termsKey]);
}

/** Clears the search; used by the chart's reset button and by "remove all filters". */
export function clearFeatureSearch(config: GeneNetworkConfig) {
    runInAction(() => {
        config.search.query = "";
        config.search.terms = [];
    });
}
