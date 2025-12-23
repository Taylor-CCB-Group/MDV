import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import type DataStore from "../../datastore/DataStore";
import { BaseReactChart } from "./BaseReactChart";
import { GeneNetworkChartComponent } from "./GeneNetworkChartComponent";
import { g, isArray } from "@/lib/utils";

export type GeneNetworkConfig = {
    mode?: "filtered" | "observableFields";
    maxGenes?: number;
    autoScroll?: boolean;
} & BaseConfig;

const DEFAULT_MAX_GENES = 100;

class GeneNetworkChartWrapper extends BaseReactChart<GeneNetworkConfig> {
    constructor(
        dataStore: DataStore,
        div: string | HTMLDivElement,
        config: GeneNetworkConfig & BaseConfig,
    ) {
        if (!config.mode) config.mode = "filtered";
        if (!config.maxGenes) config.maxGenes = DEFAULT_MAX_GENES;
        super(dataStore, div, config, GeneNetworkChartComponent);
    }
    getSettings() {
        const c = this.config;
        if (!isArray(this.config.param)) throw "expected param array";
        return [
            ...super.getSettings(),
            // g({
            //     type: "dropdown",
            //     label: "Display Mode",
            //     current_value: c.mode || "filtered",
            //     values: [
            //         [
            //             { name: "Filtered (with highlights)", value: "filtered" },
            //             { name: "ObservableFields-like", value: "observableFields" },
            //         ],
            //         "name",
            //         "value",
            //     ],
            //     func: (v) => {
            //         c.mode = v as "filtered" | "observableFields";
            //     },
            // }),
            g({
                type: "spinner",
                label: "Max genes to show",
                current_value: c.maxGenes ?? DEFAULT_MAX_GENES,
                min: 10,
                max: 1000,
                step: 10,
                func: (v) => {
                    c.maxGenes = v;
                },
            }),
        ];
    }
}

BaseChart.types["geneNetwork"] = {
    class: GeneNetworkChartWrapper,
    name: "Gene Info",
    allow_user_add: true,
    params: [
        {
            type: ["text", "text16"],
            name: "Gene ID Column (text column containing gene identifiers)",
        },
    ]
};

export default 42;

