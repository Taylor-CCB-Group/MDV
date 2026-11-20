import { g, isArray } from "@/lib/utils";
import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import type DataStore from "../../datastore/DataStore";
import { BaseReactChart } from "./BaseReactChart";
import { FEATURE_INFO_SOURCES, type FeatureInfoSource, ORGANISMS, type Organism } from "./featureInfoSources";
import { type FeatureSearch, isSearchActive } from "./featureSearch";
import { clearFeatureSearch, FEATURE_SEARCH_FILTER_KEY } from "./featureSearchHooks";
import { GeneNetworkChartComponent } from "./GeneNetworkChartComponent";

export type GeneNetworkConfig = {
    /** Whether to auto-scroll to highlighted genes coming from external interactions. */
    autoScroll?: boolean;
    /** Where card details come from; GeneNetwork (genes, human) by default. */
    source?: FeatureInfoSource;
    /** Species to look identifiers up in, for sources that support more than one. */
    organism?: Organism;
    /** Filters the whole view; see {@link FeatureSearch}. */
    search: FeatureSearch;
} & BaseConfig;

class GeneNetworkChartWrapper extends BaseReactChart<GeneNetworkConfig> {
    constructor(dataStore: DataStore, div: string | HTMLDivElement, config: GeneNetworkConfig & BaseConfig) {
        if (config.autoScroll === undefined) config.autoScroll = true;
        // set before the config is made observable, so later edits are tracked
        config.source ??= "genenetwork";
        config.organism ??= "human";
        config.search = { query: config.search?.query ?? "", terms: config.search?.terms ?? [] };
        super(dataStore, div, config, GeneNetworkChartComponent);
        const { query, terms } = this.config.search;
        if (query?.trim() || terms?.length) this.storedFilters.expect(FEATURE_SEARCH_FILTER_KEY);
        this.mobxAutorun(() => {
            this.resetButton.style.display = isSearchActive(this.config.search) ? "inline" : "none";
        });
    }
    removeFilter() {
        clearFeatureSearch(this.config);
    }
    getSettings() {
        const c = this.config;
        if (!isArray(this.config.param)) throw new Error("expected param array");
        return [
            ...super.getSettings(),
            g({
                type: "dropdown",
                label: "Information source",
                current_value: c.source ?? "genenetwork",
                values: [
                    Object.entries(FEATURE_INFO_SOURCES).map(([value, s]) => ({ name: s.label, value })),
                    "name",
                    "value",
                ],
                func: (v) => {
                    if (v === "genenetwork" || v === "uniprot") c.source = v;
                },
            }),
            g({
                type: "dropdown",
                label: "Species (UniProt)",
                current_value: c.organism ?? "human",
                values: [Object.entries(ORGANISMS).map(([value, o]) => ({ name: o.label, value })), "name", "value"],
                func: (v) => {
                    if (v === "human" || v === "mouse") c.organism = v;
                },
            }),
            g({
                type: "check",
                label: "Auto-scroll to highlighted genes",
                current_value: c.autoScroll ?? true,
                func: (x) => {
                    c.autoScroll = x;
                },
            }),
        ];
    }
}

BaseChart.types["geneNetwork"] = {
    class: GeneNetworkChartWrapper,
    name: "Gene / Protein Info",
    allow_user_add: true,
    params: [
        {
            type: ["text", "text16"],
            name: "Identifier column (gene symbols, Ensembl IDs or protein markers)",
        },
    ],
};
