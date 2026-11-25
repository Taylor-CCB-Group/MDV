import type { SimpleFeatureCollection } from "@deck.gl-community/editable-layers";
import type DataStore from "@/datastore/DataStore";
import { type CrossPcfOptions, computeCrossPCF, logCrossPCF } from "./crossPcf";
import { type DomainSource, selectionDomainSource } from "./domainSource";

/**
 * Debug triggers until there is a chart / dataflow node for cross-PCF. Not the architecture: each
 * builds a `DomainSource` and calls `computeCrossPCF`, which is what a real consumer would do.
 */
export function runCrossPcfDebug(dataStore: DataStore, source: DomainSource, opts?: CrossPcfOptions) {
    return computeCrossPCF(dataStore, source, opts)
        .then((res) => {
            logCrossPCF(res);
            return res;
        })
        .catch((e) => console.error("[cross-PCF]", e));
}

type SelectionChart = {
    dataStore: DataStore;
    config: {
        selectionFeatureCollection?: SimpleFeatureCollection;
        param: readonly unknown[];
        region?: string;
        title?: string;
    };
};

/** `__crossPcf.selection(mdv.debugChart)` in the console, after "debug / report chart" on a scatter
 *  chart with a lasso / polygon / rectangle drawn. */
function selection(chart: SelectionChart, opts?: CrossPcfOptions) {
    const source = selectionDomainSource(chart.config);
    if (!source) {
        console.warn("[cross-PCF] chart has no selection polygon");
        return;
    }
    return runCrossPcfDebug(chart.dataStore, source, opts);
}

Object.assign(globalThis, { __crossPcf: { selection } });
