import type { Param } from "@/lib/columnTypeHelpers";
import type DataStore from "../datastore/DataStore";
import type BaseChart from "./BaseChart";
import type { ExtraControl, GuiSpecType } from "./charts";


/**
 * Describes how a chart will be displayed in the 'add chart' dialog etc.
 */
export type ChartType<T extends BaseChart> = {
    /** A class extending BaseChart */
    class: new (
        dataStore: DataStore,
        contentDiv: HTMLDivElement,
        config: any,
    ) => T;
    /** The human-readable name that will appear in the 'add chart' dialog etc. */
    name: string;

    required?: string[] | ((ds: DataStore) => boolean);
    init?: (config: any, dataSource: any, extraControls: any) => void;
    extra_controls?: (dataStore: DataStore) => ExtraControl<GuiSpecType>[]; //not GuiSpec<GuiSpecType>[]; - AddChartDialog & SettingsDialog behave differently - would like to unify
    params?: { type: Param | Param[]; name: string }[];
    configEntriesUsingColumns?: string[];
    methodsUsingColumns?: string[]; // (keyof T)[]; //would like better typing here
    allow_user_add?: boolean;
};

export type ChartTypeMap = Record<string, ChartType<BaseChart>>;

/** A dictionary of all the chart types, also accessible as `BaseChart.types` */
export const chartTypes: ChartTypeMap = {};
