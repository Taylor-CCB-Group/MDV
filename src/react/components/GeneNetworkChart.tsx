import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import type DataStore from "../../datastore/DataStore";
import { BaseReactChart } from "./BaseReactChart";
import { GeneNetworkChartComponent } from "./GeneNetworkChartComponent";
import { g, isArray } from "@/lib/utils";

export type GeneNetworkConfig = {
    mode?: "filtered" | "observableFields";
    /** Whether to auto-scroll to highlighted genes coming from external interactions. */
    autoScroll?: boolean;
} & BaseConfig;

class GeneNetworkChartWrapper extends BaseReactChart<GeneNetworkConfig> {
    constructor(
        dataStore: DataStore,
        div: string | HTMLDivElement,
        config: GeneNetworkConfig & BaseConfig,
    ) {
        if (!config.mode) config.mode = "filtered";
        if (config.autoScroll === undefined) config.autoScroll = true;
        super(dataStore, div, config, GeneNetworkChartComponent);
    }
    getSettings() {
        const c = this.config;
        if (!isArray(this.config.param)) throw new Error("expected param array");
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
    name: "Gene Info",
    allow_user_add: true,
    params: [
        {
            type: ["text", "text16"],
            name: "Gene ID Column (text column containing gene identifiers)",
        },
    ]
};
