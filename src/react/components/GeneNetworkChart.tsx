import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import type DataStore from "../../datastore/DataStore";
import { BaseReactChart } from "./BaseReactChart";
import { GeneNetworkChartComponent } from "./GeneNetworkChartComponent";
import { g, isArray } from "@/lib/utils";

export type GeneNetworkConfig = {
    displayMode?: "highlighted" | "filtered";
} & BaseConfig;

class GeneNetworkChartWrapper extends BaseReactChart<GeneNetworkConfig> {
    constructor(
        dataStore: DataStore,
        div: string | HTMLDivElement,
        config: GeneNetworkConfig & BaseConfig,
    ) {
        if (!config.displayMode) config.displayMode = "highlighted";
        super(dataStore, div, config, GeneNetworkChartComponent);
    }
    getSettings() {
        const c = this.config;
        if (!isArray(this.config.param)) throw "expected param array";
        return [
            ...super.getSettings(),
            g({
                type: "dropdown",
                label: "Display Mode",
                current_value: c.displayMode || "highlighted",
                values: [
                    [
                        { name: "Highlighted", value: "highlighted" },
                        { name: "Filtered", value: "filtered" },
                    ],
                    "name",
                    "value",
                ],
                func: (v) => {
                    c.displayMode = v as "highlighted" | "filtered";
                },
            })
        ];
    }
}

BaseChart.types["geneNetwork"] = {
    class: GeneNetworkChartWrapper,
    name: "Gene Network Info",
    allow_user_add: true,
    params: [
        {
            type: ["text", "text16"],
            name: "Gene ID Column (text column containing gene identifiers)",
        },
    ]
};

export default 42;

