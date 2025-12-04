import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import { observer } from "mobx-react-lite";
import type DataStore from "@/datastore/DataStore";
import TableChartReactComponent from "./TableChartReactComponent";
import { g } from "@/lib/utils";
import { action, extendObservable, runInAction } from "mobx";

const TableChartComponent = observer(() => {
    return <TableChartReactComponent />;
});

// todo: Add params and other options to config
export type TableChartReactConfig = BaseConfig & {
    type: "table_chart_react";
    include_index?: boolean;
    column_widths?: Record<string, number>;
    sort?: { columnId: string; ascending: boolean; } | undefined;
};

// todo: Add required functions related to the features
class TableChartReact extends BaseReactChart<TableChartReactConfig> {
    constructor(
        dataStore: DataStore,
        div: HTMLDivElement,
        config: TableChartReactConfig,
    ) {
        if (!config.column_widths) config.column_widths = {};
        if (config.include_index === undefined || config.include_index === null) config.include_index = false;

        super(dataStore, div, config, TableChartComponent);

        // Extending config to make config.sort observable
        // Doing this is required as it's initialized with undefined and doesn't work if initialized like other properties
        extendObservable(this.config, {
            sort: undefined as { columnId: string; ascending: boolean } | undefined,
        });
    }

    getSettings() {
        const settings = super.getSettings();
        const c = this.config;
        return [
            ...settings,
            g({
                type: "check",
                current_value: c.include_index || false,
                label: "Include Index Column",
                func: action((x: boolean) => {
                    c.include_index = x;
                }),
            }),
        ];
    }
}

// todo: Sync the options of config including the params
BaseChart.types["table_chart_react"] = {
    name: "Table Chart (React)",
    class: TableChartReact,
    allow_user_add: true,
    params: [
        {
            type: "_multi_column:all",
            name: "Columns To Display",
        },
    ],
    extra_controls: () => [
        {
            type: "check",
            name: "include_index",
            label: "Include Index",
        },
    ],
    init: (config, ds, ec) => {
        config.include_index = ec.include_index;
    },
};

// export for side effect of HMR
export default 42;


