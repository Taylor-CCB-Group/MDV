import BaseChart from "../../charts/BaseChart";
import { BaseReactChart, type BaseConfig } from "./BaseReactChart";
import SelectionDialogComponent from "./SelectionDialogComponent";

export type CategoryFilter = { category: string[], operand?: "or" | "and" };

export type SelectionDialogFilter = CategoryFilter | [number, number];

export type SelectionDialogConfig = {
    type: "selection_dialog_experimental";
    filters?: Record<string, SelectionDialogFilter | null>;
};

class SelectionDialogReact extends BaseReactChart<SelectionDialogConfig> {
    constructor(dataStore, div, config: SelectionDialogConfig & BaseConfig) {
        if (!config.filters) {
            config.filters = {};
            for (const col of config.param) {
                config.filters[col] = null;
            }
        }
        for (const col of config.param) {
            if (!config.filters[col]) {
                config.filters[col] = null;
            }
        }
        // makeAutoObservable(config); //super will do this 
        super(dataStore, div, config, SelectionDialogComponent);
        console.log("todo add reset button & associated logic in SelectionDialogReact");
    }
}

BaseChart.types["selection_dialog_experimental"] = {
    name: "Selection Dialog (Experimental)",
    class: SelectionDialogReact,
    params: [
        {
            type: "_multi_column:all",
            name: "Columns To filter",
        },
    ],
}