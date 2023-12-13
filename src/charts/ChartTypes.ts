import DataStore from "../datastore/DataStore";
import BaseChart from "./BaseChart";
import type { DataSource, GuiSpec } from "./charts";

export type ChartType<T extends BaseChart> = {
    /** A class extending BaseChart */
    "class": new (dataStore: DataStore, contentDiv: HTMLDivElement, config: any) => T;
    /** The human-readable name that will appear in the 'add chart' dialog etc. */
    name: string;
    
    required?: string[] | ((ds: DataSource) => unknown); //TODO: better typing here (& there & everywhere)
    init?: (config: any, dataSource: any, extraControls: any) => void;
    extra_controls?: (dataSource: DataSource) => GuiSpec<any>[];
    params?: { type: string | string[], name: string }[];
    configEntriesUsingColumns?: string[];
    methodsUsingColumns?: string[]; // (keyof T)[]; //would like better typing here
}

export type ChartTypeMap = Record<string, ChartType<BaseChart>>;

/** A dictionary of all the chart types, also accessible as `BaseChart.types` */
export const chartTypes: ChartTypeMap = {};