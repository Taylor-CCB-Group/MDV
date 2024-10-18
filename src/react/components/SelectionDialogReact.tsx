import { action } from "mobx";
import BaseChart from "../../charts/BaseChart";
import { BaseReactChart, type BaseConfig } from "./BaseReactChart";
import SelectionDialogComponent from "./SelectionDialogComponent";
import { observer } from "mobx-react-lite";

export type CategoryFilter = { category: string[] };
export type MultiTextFilter = CategoryFilter & { operand: "or" | "and" };
export type RangeFilter = [number, number];
export type UniqueFilter = string;
export type SelectionDialogFilter = CategoryFilter | MultiTextFilter | UniqueFilter | RangeFilter;

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
        //nb, considered `this.mobxAutorun` for showing/hiding reset button, but we use a hook.
        super(dataStore, div, config, observer(SelectionDialogComponent));
    }
    removeFilter(): void {
        action(() => {
            for (const key in this.config.filters) {
                this.config.filters[key] = null;
            }
        })();
    }
}

BaseChart.types["selection_dialog"] = {
    name: "Selection Dialog",
    class: SelectionDialogReact,
    params: [
        {
            type: "_multi_column:all",
            name: "Columns To filter",
        },
    ],
}
// BaseChart.types["selection_dialog_experimental"] = BaseChart.types["selection_dialog"];