import { action } from "mobx";
import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import type DataStore from "../../datastore/DataStore";
import { BaseReactChart } from "./BaseReactChart";
import { HighlightedFeatureComponent } from "./HighlightedFeatureComponent";
import type { GuiSpec } from "@/charts/charts";

type HighlightedFeatureConfig = {
    text: string;
};

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
            {
                type: "dropdown",
                // name: "Column",
                label: "Column",
                current_value: this.config.param[0],
                values: [this.dataSource.dataStore.columns.map((c) => c.name)],
                func: action((v: string) => {
                    this.config.param[0] = v;
                }),
            },
            {
                label: "Markdown Text",
                // name: "text",
                type: "textbox",
                current_value: c.text || "",
                func: action((x) => {
                    c.text = x;
                }),
            },
        ] as GuiSpec<any>[];
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
