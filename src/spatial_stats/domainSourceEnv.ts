import type BaseChart from "@/charts/BaseChart";
import type { BaseConfig } from "@/charts/BaseChart";
import type DataStore from "@/datastore/DataStore";
import { getOrCreateGateManager } from "@/react/gates/useGateManager";
import type { DomainSource } from "./domainSource";
import type { DomainSourceEnv } from "./domainSourceRef";
import { regionDomainSource } from "./regionDomain";
import { spatialRegionKeys } from "./regionSpatialData";

/** Stable whole-region sources per datastore, so consumers can key work on their identity. */
const regionSources = new WeakMap<DataStore, Map<string, DomainSource>>();
function regionSourceFor(dataStore: DataStore, region?: string) {
    let byKey = regionSources.get(dataStore);
    if (!byKey) regionSources.set(dataStore, (byKey = new Map()));
    const key = region ?? "";
    let s = byKey.get(key);
    if (!s) byKey.set(key, (s = regionDomainSource(dataStore, region)));
    return s;
}

/** The live app's sources for `dataStore`: every chart that might host a selection, and its gates. */
export function domainSourceEnv(dataStore: DataStore): DomainSourceEnv {
    const charts = Object.values(window.mdv.chartManager.charts).map((c: { chart: BaseChart<BaseConfig> }) => ({
        dataStore: c.chart.dataStore,
        config: c.chart.config,
    }));
    return {
        dataStore,
        charts,
        gates: getOrCreateGateManager(dataStore).gates,
        regions: spatialRegionKeys(dataStore),
        regionSource: (region) => regionSourceFor(dataStore, region),
    };
}
