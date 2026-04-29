import BaseChart, { type BaseConfig } from "@/charts/BaseChart";
import type DataStore from "@/datastore/DataStore";
import { g } from "@/lib/utils";
import { BaseReactChart } from "./BaseReactChart";
import TextBoxChartComponent from "./textbox/TextBoxChartComponent";

export type TextBoxChartConfig = BaseConfig & {
    type: "text_box_chart" | "text_box";
    text: string;
};

class TextBoxChartReactWrapper extends BaseReactChart<TextBoxChartConfig> {
    constructor(
        dataStore: DataStore,
        div: string | HTMLDivElement,
        config: TextBoxChartConfig,
    ) {
        if (!config.text) config.text = "";
        super(dataStore, div, config, TextBoxChartComponent);
    }

    getSettings() {
        const c = this.config;
        return [
            ...super.getSettings(),
            g({
                label: "Markdown Text",
                type: "textbox",
                current_value: c.text ?? "",
                func: (text) => {
                    c.text = text;
                },
            }),
        ];
    }
}

BaseChart.types["text_box_chart"] = {
    class: TextBoxChartReactWrapper,
    name: "Text Box",
    params: [],
    extra_controls: () => [
        {
            type: "textbox",
            name: "text",
            label: "Markdown Text:",
        },
    ],
    init: (config, _dataSource, extraControls) => {
        config.text = extraControls.text ?? "";
    },
};
