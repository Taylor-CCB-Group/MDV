import { action } from "mobx";
import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import SelectionDialogComponent from "./SelectionDialogComponent";
import type DataStore from "@/datastore/DataStore";
import { flattenFields } from "@/lib/columnTypeHelpers";
import { g } from "@/lib/utils";
import type { ChartColumnImpact } from "@/types/columnRemovalTypes";

// export type CommonFilter = { noClear?: boolean, invert?: boolean }; //todo - implement these features.
export type CategoryFilter = { category: string[] };
export type MultiTextFilter = CategoryFilter & { operand: "or" | "and" };
export type RangeFilter = [number, number];
export type UniqueFilter = string;
export type SelectionDialogFilter = (CategoryFilter | MultiTextFilter | UniqueFilter | RangeFilter);// & CommonFilter;

export type SelectionDialogConfig = {
    type: "selection_dialog";
    noClearFilters?:boolean;
    filters: Record<string, SelectionDialogFilter | null>;
    order?: Record<string, number>;
} & BaseConfig;

class SelectionDialogReact extends BaseReactChart<SelectionDialogConfig> {
    constructor(dataStore: DataStore, div: HTMLDivElement, config: SelectionDialogConfig & BaseConfig) {
        if (!config.filters) {
            config.filters = {};
            for (const col of config.param) {
                //@ts-expect-error MultiColumnQuery cannot be used as index
                config.filters[col] = null;
            }
        }

        if (!config.order) {
            config.order = {};
        }
        //for legacy configs, if left undefined makeAutoObservable will 
        // not work correctly.
        if (config.noClearFilters === undefined) {
            config.noClearFilters = false;
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

    onColumnRemoved(column: string, impact?: ChartColumnImpact): boolean {
        if (impact?.action === "delete") {
            return super.onColumnRemoved(column, impact);
        }
        const currentParam = Array.isArray(this.config.param) ? this.config.param : [];
        const nextParam =
            impact?.nextParam ?? currentParam.filter((field) => !flattenFields(field).includes(column));
        // `config.param` can contain query-backed FieldSpecs, so diff the
        // flattened concrete field names rather than assuming `column` is the
        // only key that disappears from filters/order.
        const nextFieldSet = new Set(nextParam.flatMap((field) => flattenFields(field)));
        const removedFields = currentParam
            .flatMap((field) => flattenFields(field))
            .filter((field) => !nextFieldSet.has(field));

        const chartImpact: ChartColumnImpact = impact ?? {
            // Fallback for any non-standard caller. In the normal flow
            // ChartManager provides a richer impact object.
            chartId: this.config.id,
            chartTitle: this.config.title,
            chartType: this.config.type,
            chartTypeLabel: BaseChart.types[this.config.type]?.name ?? this.config.type,
            isSourceChart: false,
            action: "update",
            reasons: [],
            paramSlotImpacts: [],
            nextParam,
            clearColorBy: false,
            clearBackgroundFilter: false,
        };

        // Delegate param pruning and shared config cleanup to BaseChart instead
        // of mutating `config.param` directly here.
        const didDelete = super.onColumnRemoved(column, {
            ...chartImpact,
            action: "update",
            nextParam,
        });
        if (didDelete) {
            return true;
        }
        if (removedFields.length === 0) {
            // BaseChart already handled the generic cleanup; there are just no
            // selection-dialog-specific filter/order keys to prune.
            return false;
        }

        action(() => {
            // The shared analyzer has already decided that this chart survives;
            // this override just keeps filter/order state aligned with pruned params.
            for (const field of removedFields) {
                delete this.config.filters[field];
            }

            const nextOrder = { ...(this.config.order ?? {}) };
            // `order` is stored as dense integer positions, so after removing
            // one or more fields every later entry must shift left.
            const removedOrders = removedFields
                .map((field) => nextOrder[field])
                .filter((order): order is number => order !== undefined)
                .sort((a, b) => a - b);
            for (const field of removedFields) {
                delete nextOrder[field];
            }
            if (removedOrders.length > 0) {
                for (const field in nextOrder) {
                    const currentOrder = nextOrder[field];
                    const shift = removedOrders.filter((order) => order < currentOrder).length;
                    if (shift > 0) {
                        nextOrder[field] -= shift;
                    }
                }
            }
            this.config.order = nextOrder;
        })();

        return false;
    }

    getSettings(){
          const settings = super.getSettings();
          const c = this.config;
          return [
            ...settings,
            g({
                //set whether filters are cleared on Reset All
                type: "check",
                current_value: c.noClearFilters || false,
                label: "Filters remain on Reset All",
                func: (x) => {
                    // func() is already called in an action, so no need to wrap again
                    c.noClearFilters = x;
                }
            })
          ]
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

export default SelectionDialogReact;
