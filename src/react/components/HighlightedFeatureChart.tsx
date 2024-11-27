import { action } from "mobx";
import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import type DataStore from "../../datastore/DataStore";
import { BaseReactChart } from "./BaseReactChart";
import { HighlightedFeatureComponent } from "./HighlightedFeatureComponent";
import { g } from "@/lib/utils";

type HighlightedFeatureConfig = {
    text: string;
} & BaseConfig;

class HighlightedFeatureChartWrapper extends BaseReactChart<HighlightedFeatureConfig> {
    constructor(
        dataStore: DataStore,
        div: string | HTMLDivElement,
        config: HighlightedFeatureConfig & BaseConfig,
    ) {
        if (!config.text) config.text = "";
        super(dataStore, div, config, HighlightedFeatureComponent);
    }
    getSettings() {
        const c = this.config;
        return [
            ...super.getSettings(),
            g({
                type: "dropdown",
                // name: "Column",
                label: "Column",
                //@ts-expect-error needs fixing
                current_value: this.config.param[0],
                values: [this.dataSource.dataStore.columns.map((c) => c.name)],
                func: action((v: string) => {
                    //@ts-expect-error needs fixing
                    this.config.param[0] = v;
                }),
            }),
            g({
                label: "Markdown Text",
                // name: "text",
                type: "textbox",
                current_value: c.text || "",
                func(x) {
                    c.text = x;
                },
            }),
        ];
    }
}

BaseChart.types["highlightedFeature"] = {
    class: HighlightedFeatureChartWrapper,
    name: "Highlighted Feature",
    allow_user_add: false,
    // required: ['muspan'],
    params: [
        {
            type: "_multi_column:all",
            name: "Columns to include in server query",
        },
    ],
    extra_controls: () => [
        {
            type: "textbox",
            name: "text",
            label: "Markdown Text:",
        },
    ],
    init: (config, ds, extraControls) => {
        config.text = extraControls.text || "";
    },
};

export default 42;
