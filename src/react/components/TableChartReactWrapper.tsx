import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import { observer } from "mobx-react-lite";
import type DataStore from "@/datastore/DataStore";
import TableChartReactComponent from "./TableChartReactComponent";

const TableChartComponent = observer(() => {
    return <TableChartReactComponent />;
});

// todo: Add params and other options to config
export type TableChartReactConfig = BaseConfig & {
    type: "table_chart_react";
};

// todo: Add required functions related to the features
class TableChartReact extends BaseReactChart<TableChartReactConfig> {
    constructor(
        dataStore: DataStore,
        div: HTMLDivElement,
        config: TableChartReactConfig,
    ) {
        super(dataStore, div, config, TableChartComponent);
    }
}

// todo: Sync the options of config including the params
BaseChart.types["table_chart_react"] = {
    name: "Table Chart (React)",
    class: TableChartReact,
    allow_user_add: true,
    params: [],
};

// export for side effect of HMR
export default 42;


