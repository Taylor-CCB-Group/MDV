import type DataStore from "../datastore/DataStore";
import type BaseChart from "./BaseChart";
import type { BaseConfig } from "./BaseChart";
import type { DataType, ExtraControl, GuiSpecType } from "./charts";

//new Set(Object.values(BaseChart.types).flatMap(t => t.params).filter(Boolean).flatMap(p => p.type))
/** annotation of what kind of column(s) a given param will accept */
// export type Param = "text" | "number" | "multitext" | "text16" | "_multi_column:number" | "_multi_column:all";
export type Param = "text" | "number" | "multitext" | "text16" | "_multi_column:number" | "_multi_column:all";

//chatGPT to the rescue
export type BaseChartConstructor<TProps extends BaseConfig> = new (...args: any[]) => BaseChart<TProps>;
/**
 * Describes how a chart will be displayed in the 'add chart' dialog etc.
 * 
 * This is something we may want to turn into a more formal JSON schema / zod validated etc.
 */
export type ChartType<T extends BaseConfig> = {
    /** A class extending BaseChart */
    class: BaseChartConstructor<T>;
    /** The human-readable name that will appear in the 'add chart' dialog etc. */
    name: string;

    /** 
     * If particular features need to be enabled on the `DataStore`, these can be listed here, or a function can be provided
     * to indicate if the given `DataStore` is suitable.
     */
    required?: string[] | ((ds: DataStore) => boolean);
    init?: (config: any, dataSource: any, extraControls: any) => void;
    extra_controls?: (dataStore: DataStore) => ExtraControl<GuiSpecType>[]; //not GuiSpec<GuiSpecType>[]; - AddChartDialog & SettingsDialog behave differently - would like to unify
    params?: { type: Param | Param[]; name: string; reactive?: boolean }[];
    configEntriesUsingColumns?: string[];
    /**
     * @deprecated - use {@link loadColumnData} decorator directly on the method.
     */
    methodsUsingColumns?: string[]; // (keyof T)[]; //would like better typing here
    allow_user_add?: boolean;
    // given that any chart with color_by has already been made to work with column queries, 
    // we probably want to flag this on a more granular level...
    // Things like 'color_by' aren't represented here anyway (and can't be set in the AddChartDialog as of now)
    // allowDynamicColumns?: boolean;
};

export type ChartTypeMap = Record<string, ChartType<any>>;

/** A dictionary of all the chart types, also accessible as `BaseChart.types` */
export const chartTypes: ChartTypeMap = {};
