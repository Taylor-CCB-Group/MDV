import { g } from "@/lib/utils";
import BaseChart, { type BaseConfig } from "../charts/BaseChart";
import type DataStore from "../datastore/DataStore";
import { BaseReactChart } from "../react/components/BaseReactChart";
import { CrossPcfChartComponent } from "./CrossPcfChartComponent";
import { type LabelPair, togglePair } from "./CrossPcfMatrix";
import { domainSourceEnv } from "./domainSourceEnv";
import { type DomainSourceRef, listDomainSources, refKey } from "./domainSourceRef";

export type CrossPcfConfig = BaseConfig & {
    /** Categorical column whose categories are the cell types. */
    type_column: string | null;
    /** Pairs drawn as g(r) curves and outlined in the matrix, as category labels of `type_column`:
     *  anchor A → neighbour B (g is asymmetric). Labels, not codes, so they survive a domain with fewer types. */
    pairs: LabelPair[];
    /** What the observation window is; unset → a drawn selection, else the whole region. */
    domain_source: DomainSourceRef | null;
    /** In the analysis frame's unit. */
    r_max: number;
    n_bins: number;
    edge_corrected: boolean;
    /** Bin shown in the matrix, or the mean over r. */
    r_bin: number | "mean";
    /** µm per unit, applied only when the store leaves its unit unstated (0 = none). */
    um_per_unit: number;
};

class CrossPcfChart extends BaseReactChart<CrossPcfConfig> {
    constructor(dataStore: DataStore, div: string | HTMLDivElement, config: CrossPcfConfig) {
        config.r_max ??= 100;
        config.n_bins ??= 20;
        config.edge_corrected ??= true;
        config.r_bin ??= "mean";
        config.um_per_unit ??= 0;
        // Unset is `null`, not absent: the config is JSON-copied then made observable key by key, so a key
        // that isn't there (or is `undefined`) when that happens is never tracked.
        config.domain_source ??= null;
        config.type_column ??= null;
        config.pairs ??= [];
        // Semantic props only; BaseChart still expects the array.
        config.param ??= [];
        super(dataStore, div, config, CrossPcfChartComponent);
    }
    getSettings() {
        const c = this.config;
        const sources = listDomainSources(domainSourceEnv(this.dataStore));
        const values = [{ name: "(a drawn selection, else the whole region)", value: "" }].concat(
            sources.map((s) => ({ name: s.label, value: refKey(s.ref) })),
        );
        const categories = c.type_column ? (this.dataStore.columnIndex[c.type_column]?.values ?? []) : [];
        // The A/B pickers only stage a pair; the button commits it to `pairs`.
        const last = c.pairs.at(-1);
        const pending = { a: last?.a ?? "", b: last?.b ?? "" };
        const categoryValues = [{ name: "(none)", value: "" }].concat(
            categories.map((v: string) => ({ name: v, value: v })),
        );
        return [
            ...super.getSettings(),
            g({
                type: "dropdown",
                label: "Domain",
                current_value: c.domain_source ? refKey(c.domain_source) : "",
                values: [values, "name", "value"],
                func: (v) => {
                    c.domain_source = sources.find((s) => refKey(s.ref) === v)?.ref ?? null;
                },
            }),
            g({
                type: "column",
                label: "Cell type column",
                current_value: c.type_column ?? "",
                columnType: ["text", "text16"],
                func: (x) => {
                    if (typeof x !== "string" || x === c.type_column) return;
                    c.type_column = x;
                    c.pairs = [];
                },
            }),
            g({
                type: "dropdown",
                label: "Category A (anchor)",
                current_value: pending.a,
                values: [categoryValues, "name", "value"],
                func: (v) => {
                    pending.a = v;
                },
            }),
            g({
                type: "dropdown",
                label: "Category B (neighbour)",
                current_value: pending.b,
                values: [categoryValues, "name", "value"],
                func: (v) => {
                    pending.b = v;
                },
            }),
            g({
                type: "button",
                label: "Add / remove pair A → B",
                current_value: null,
                func: () => {
                    if (pending.a && pending.b) c.pairs = togglePair(c.pairs, pending.a, pending.b);
                },
            }),
            g({
                type: "button",
                label: `Clear pairs (${c.pairs.length})`,
                current_value: null,
                func: () => {
                    c.pairs = [];
                },
            }),
            g({
                type: "spinner",
                label: "r max",
                current_value: c.r_max,
                min: 1,
                max: 100_000,
                step: 1,
                func: (v) => {
                    c.r_max = v;
                },
            }),
            g({
                type: "spinner",
                label: "Bins",
                current_value: c.n_bins,
                min: 1,
                max: 200,
                step: 1,
                func: (v) => {
                    c.n_bins = v;
                    if (typeof c.r_bin === "number" && c.r_bin >= v) c.r_bin = "mean";
                },
            }),
            g({
                type: "check",
                label: "Edge correction (clip annuli to the domain)",
                current_value: c.edge_corrected,
                func: (x) => {
                    c.edge_corrected = x;
                },
            }),
            g({
                type: "spinner",
                label: "µm per unit, if the store doesn't say (0 = none)",
                current_value: c.um_per_unit,
                min: 0,
                max: 1e6,
                step: 0.0001,
                func: (v) => {
                    c.um_per_unit = v;
                },
            }),
        ];
    }
}

BaseChart.types["cross_pcf"] = {
    class: CrossPcfChart,
    name: "Cross-PCF matrix (prototype)",
    allow_user_add: true,
    params: [],
    extra_controls: (dataStore: DataStore) => {
        const values = dataStore
            .getColumnList(["text", "text16"])
            .map((c) => ({ name: c.name ?? "", value: c.field ?? "" }));
        const preferred = dataStore.regions?.default_color;
        return [
            {
                type: "dropdown",
                name: "type_column",
                label: "Cell type column",
                values,
                defaultVal: values.some((v) => v.value === preferred) ? preferred : values[0]?.value,
            },
        ];
    },
};
