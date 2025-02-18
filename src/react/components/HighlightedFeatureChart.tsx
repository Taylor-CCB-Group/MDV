import { action } from "mobx";
import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import type DataStore from "../../datastore/DataStore";
import { BaseReactChart } from "./BaseReactChart";
import { HighlightedFeatureComponent } from "./HighlightedFeatureComponent";
import { g, isArray } from "@/lib/utils";

export type HighlightedFeatureConfig = {
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
        if (!isArray(this.config.param)) throw "expected param array";
        return [
            ...super.getSettings(),
            g({
                // untested / chart is generally unused
                type: "multicolumn",
                label: "Column",
                current_value: this.config.param,
                func: (v) => {
                    this.config.param = v;
                },
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
