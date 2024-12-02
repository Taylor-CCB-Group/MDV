import { action } from "mobx";
import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import SelectionDialogComponent from "./SelectionDialogComponent";
import { observer } from "mobx-react-lite";
import type DataStore from "@/datastore/DataStore";
import { g } from "@/lib/utils";

export type CategoryFilter = { category: string[] };
export type MultiTextFilter = CategoryFilter & { operand: "or" | "and" };
export type RangeFilter = [number, number];
export type UniqueFilter = string;
export type SelectionDialogFilter = CategoryFilter | MultiTextFilter | UniqueFilter | RangeFilter;

export type SelectionDialogConfig = {
    type: "selection_dialog";
    filters: Record<string, SelectionDialogFilter | null>;
} & BaseConfig;

class SelectionDialogReact extends BaseReactChart<SelectionDialogConfig> {
    constructor(dataStore: DataStore, div: HTMLDivElement, config: SelectionDialogConfig & BaseConfig) {
        if (!config.filters) {
            config.filters = {};
            //@ts -ignore ! @ts-expect-error is inconsistent between editor & cli???
            for (const col of config.param) {
                //@ts-expect-error MultiColumnQuery cannot be used as index
                config.filters[col] = null;
            }
        }
        for (const col of config.param) {
            //@ts-expect-error MultiColumnQuery cannot be used as index
            if (!config.filters[col]) {
                //@ts-expect-error MultiColumnQuery cannot be used as index
                config.filters[col] = null;
            }
        }
        // makeAutoObservable(config); //super will do this
        //nb, considered `this.mobxAutorun` for showing/hiding reset button, but we use a hook.
        super(dataStore, div, config, SelectionDialogComponent);
    }
    removeFilter(): void {
        action(() => {
            for (const key in this.config.filters) {
                this.config.filters[key] = null;
            }
        })();
    }
    getSettings() {
        // todo: add settings widget for 'column' with some properties somewhat similar to params type.
        const settings = super.getSettings();
        // >>> it is indeed broken at runtime - need general development of multi-column settings
        // and how they relate to the mobx config store.
        settings.push(g({
            //!!this should be properly implemented...
            type: "multicolumn",
            label: "Columns To filter",
            current_value: this.config.param,
            func: (v) => {
                this.config.param = v;
                // this.removeFilter();
            }
        }));
        return settings;
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
