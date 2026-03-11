import BaseChart, { type BaseConfig } from "../../charts/BaseChart";
import { BaseReactChart } from "./BaseReactChart";
import type DataStore from "@/datastore/DataStore";
import TableChartReactComponent from "./TableChartReactComponent";
import { g } from "@/lib/utils";
import { action } from "mobx";
import type { SlickgridReactInstance } from "slickgrid-react";
import { getTableExportBlob } from "@/datastore/dataExportUtils";
import { createEl } from "@/utilities/ElementsTyped";

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

class TableChartReact extends BaseReactChart<TableChartReactConfig> {
    gridRef?: { current: SlickgridReactInstance | null };
    private downloadIconSpan: HTMLSpanElement | null = null;

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
            const blob = await getTableExportBlob(this.dataStore, columns, {
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

    getConfig() {
        const config = super.getConfig();
        
        // Serialize current grid widths to config.column_widths for persistence
        // This ensures persisted state reflects actual grid state at serialization time
        // Grid manages widths internally during runtime, but we serialize them here for persistence
        if (this.gridRef?.current?.slickGrid) {
            try {
                const columns = this.gridRef.current.slickGrid.getColumns();
                const columnWidths: Record<string, number> = {};
                
                for (const col of columns) {
                    if (col.width && col.width !== 100) {
                        columnWidths[col.field] = col.width;
                    }
                }
                
                config.column_widths = columnWidths;
            } catch (e) {
                // If grid is not available or error occurs, keep existing config.column_widths
                console.warn("Could not read grid widths in getConfig():", e);
            }
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
        // nb - may be undefined if user hasn't touched that control
        config.include_index = ec.include_index;
    },
};

// export for side effect of HMR
export default 42;
