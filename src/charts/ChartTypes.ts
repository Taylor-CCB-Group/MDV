import type DataStore from "../datastore/DataStore";
import type BaseChart from "./BaseChart";
import type { BaseConfig } from "./BaseChart";
import type { DataType, ExtraControl, GuiSpecType } from "./charts";

//new Set(Object.values(BaseChart.types).flatMap(t => t.params).filter(Boolean).flatMap(p => p.type))
/** annotation of what kind of column type a given param will accept */
// export type Param = "text" | "number" | "multitext" | "text16" | "_multi_column:number" | "_multi_column:all";
export type Param = "text" | "number" | "multitext" | "text16" | "_multi_column:number" | "_multi_column:all";

//chatGPT to the rescue
type BaseChartConstructor<TProps extends BaseConfig> = new (...args: any[]) => BaseChart<TProps>;
/**
 * Describes how a chart will be displayed in the 'add chart' dialog etc.
 */
export type ChartType<T extends BaseConfig> = {
    /** A class extending BaseChart */
    class: BaseChartConstructor<T>;
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

export type ChartTypeMap = Record<string, ChartType<any>>;

/** A dictionary of all the chart types, also accessible as `BaseChart.types` */
export const chartTypes: ChartTypeMap = {};
