import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import type DataStore from "@/datastore/DataStore";
import TableChartReactComponent from "./TableChartReactComponent";
import { g } from "@/lib/utils";
import { action } from "mobx";
import type { SlickgridReactInstance } from "slickgrid-react";
import { getTableExportBlob } from "@/datastore/dataExportUtils";
import { createEl } from "@/utilities/ElementsTyped";
import { DataModel } from "@/table/DataModel";
import { flattenFields, type FieldSpec } from "@/lib/columnTypeHelpers";

const TableChartComponent = () => {
    return <TableChartReactComponent />;
};

export type TableChartReactConfig = BaseConfig & {
    type: "table_chart_react";
    include_index?: boolean;
    column_widths?: Record<string, number>;
    order?: Record<string, number>;
    sort?: { columnId: string; ascending: boolean } | null;
};

/**
 * Adapts and validates the config before it is made observable by makeAutoObservable.
 * 
 * This ensures all properties exist before makeAutoObservable is called, which is required
 * for MobX to properly track them. All properties should be initialized to their default
 * values here.
 * 
 * In future, this may be replaced with a more formal validation step (e.g., zod schema).
 * Probably incorporated into the `initialiseChartConfig` function, which will be able to lookup
 * from schemas...
 */
function adaptConfig(config: TableChartReactConfig): TableChartReactConfig {
    // Ensure all properties exist before makeAutoObservable
    // This ensures MobX can properly track them
    if (!config.column_widths) config.column_widths = {};
    if (!config.order) config.order = {};
    if (config.include_index === undefined || config.include_index === null) {
        config.include_index = false;
    }
    if (!('sort' in config)) {
        // Assigning this to null instead of undefined to make it observable
        config.sort = null;
    }
    return config;
}

export class TableChartReact extends BaseReactChart<TableChartReactConfig> {
    private gridRef?: { current: SlickgridReactInstance | null };
    private addColumnDialogOpener?: () => void;
    private downloadIconSpan: HTMLSpanElement | null = null;
    private addColumnIconSpan: HTMLSpanElement | null = null;

    constructor(dataStore: DataStore, div: HTMLDivElement, config: TableChartReactConfig) {
        // Adapt config before calling super constructor, which will make the config observable
        // This ensures all properties exist and are properly initialized before makeAutoObservable runs
        config = adaptConfig(config);
        super(dataStore, div, config, TableChartComponent);

        // Download data
        this.downloadIconSpan = this.addMenuIcon("fas fa-download", "Download data", {
            func: () => {
                this.downloadData();
            },
        }) as HTMLSpanElement;

        const isEditMode = window.mdv?.chartManager?.config?.permission === "edit";

        if (isEditMode) {
            this.addColumnIconSpan = this.addMenuIcon("fas fa-plus", "Add Column", {
                func: () => {
                    this.openAddColumnDialog();
                },
            }) as HTMLSpanElement;
        }
    }

    setGridRef(gridRef: { current: SlickgridReactInstance | null }) {
        this.gridRef = gridRef;
    }

    setAddColumnDialogOpener(opener?: () => void) {
        this.addColumnDialogOpener = opener;
    }

    openAddColumnDialog() {
        this.addColumnDialogOpener?.();
    }

    private setDownloadIcon(downloading: boolean) {
        if (!this.downloadIconSpan) return;
        const icon = this.downloadIconSpan.firstElementChild as HTMLElement | null;
        if (!icon) return;
        // Show spinner when downloading
        if (downloading) {
            icon.classList.remove("fa-download");
            icon.classList.add("fa-spinner", "fa-spin");
            this.downloadIconSpan.style.pointerEvents = "none";
        } else {
            icon.classList.remove("fa-spinner", "fa-spin");
            icon.classList.add("fa-download");
            this.downloadIconSpan.style.pointerEvents = "";
        }
    }

    async downloadData() {
        const param = this.config.param;
        // If no columns are displayed, prompt an alert
        if (!param || !Array.isArray(param) || param.length === 0) {
            alert("Add columns to the table before downloading");
            return;
        }

        // Flatten the params to columns
        const columns = param.flatMap((p) => (typeof p === "string" ? p : p.fields));
        if (columns.length === 0) {
            alert("Add columns to the table before downloading");
            return;
        }

        // Set the icon to spinner
        this.setDownloadIcon(true);
        try {
            const dataModel = new DataModel(this.dataStore, { autoupdate: false });
            dataModel.setColumns(columns);
            const blob = await getTableExportBlob(dataModel, {
                includeIndex: this.config.include_index ?? true,
            });

            const url = window.URL.createObjectURL(blob);

            const save = createEl(
                "a",
                {
                    download: `${this.dataStore.name}.txt`,
                    target: "_blank",
                    href: url,
                },
                document.body,
            );
            save.click();
            save.remove();
            setTimeout(() => window.URL.revokeObjectURL(url), 0);
        } catch (err) {
            console.error("Failed to download table data: ", err);
            alert("Failed to download table data");
        } finally {
            this.setDownloadIcon(false);
        }
    }

    // Overriding the onColumnRemoved to update table config and avoid removing the whole table
    onColumnRemoved(column: string) {
        if (!Array.isArray(this.config.param)) {
            return false;
        }
        // `config.param` can contain either plain field names or query-backed FieldSpec objects
        // flatten each spec to concrete field names so removals work for both
        const containsRemovedColumn = (fieldSpec: FieldSpec) =>
            flattenFields(fieldSpec).includes(column);
        if (!this.config.param.some(containsRemovedColumn)) {
            return false;
        }

        action(() => {
            // Remove column from config.param, config.order
            this.config.param = this.config.param.filter(
                (fieldSpec) => !containsRemovedColumn(fieldSpec),
            );
            const nextOrder = { ...(this.config.order ?? {}) };
            const removedOrder = nextOrder[column];
            delete nextOrder[column];

            // Shift the index of the columns back after the removed column by 1
            if (removedOrder !== undefined) {
                for (const field in nextOrder) {
                    if (nextOrder[field] > removedOrder) {
                        nextOrder[field] -= 1;
                    }
                }
            }
            this.config.order = nextOrder;

            // Reset sort if sorted by removed column
            if (this.config.sort?.columnId === column) {
                this.config.sort = null;
            }

            // Remove the column width of removed column
            if (this.config.column_widths?.[column] != null) {
                delete this.config.column_widths[column];
            }
        })();

        return false;
    }

    getConfig() {
        const config = super.getConfig();

        const gridColumns = this.gridRef?.current?.slickGrid?.getColumns();
        const visibleFields = gridColumns
            ?.map((column) => column.field)
            .filter(
                (field): field is string =>
                    typeof field === "string" &&
                    field !== "__index__" &&
                    Boolean(this.dataStore.columnIndex[field]),
            );

        const fallbackFields =
            Array.isArray(this.config.param)
                ? this.config.param.flatMap((fieldSpec) => flattenFields(fieldSpec))
                    .filter((field) => Boolean(this.dataStore.columnIndex[field]))
                : [];

        const persistedFields = visibleFields ?? fallbackFields;
        config.param = persistedFields;

        config.order = Object.fromEntries(
            persistedFields.map((field, index) => [field, index]),
        );

        const columnWidths: Record<string, number> = {};
        for (const column of gridColumns ?? []) {
            if (
                typeof column.field === "string" &&
                column.field !== "__index__" &&
                persistedFields.includes(column.field) &&
                column.width &&
                column.width !== 100
            ) {
                columnWidths[column.field] = column.width;
            }
        }
        config.column_widths = columnWidths;

        if (config.sort?.columnId && !persistedFields.includes(config.sort.columnId)) {
            config.sort = null;
        }

        return config;
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

BaseChart.types["table_chart_react"] = {
    name: "Table",
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
        // nb - may be undefined if user hasn't touched that control
        config.include_index = ec.include_index;
    },
};

// export for side effect of HMR
export default 42;
