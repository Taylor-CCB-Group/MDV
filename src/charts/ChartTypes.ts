import type DataStore from "../datastore/DataStore";
import type BaseChart from "./BaseChart";
import type { BaseConfig } from "./BaseChart";
import type { DataType, ExtraControl, GuiSpecType } from "./charts";

/**
 * Annotation of what kind of column (or columns) type a given param will accept.
 * 
 * Note that these are not the same as the `DataType` values used in `DataSource`. As well as allowing
 * for `"_multi_column:number"` and `"_multi_column:all"` types, we also have `"number"` - loosely
 * equivalent to `export type NumberDataType = "integer" | "double" | "int32";` in `charts.d.ts`.
 * 
 * The list of legal values is deemed exhaustive based on the evaluation of
 * ```ts
 * new Set(Object.values(BaseChart.types).flatMap(t => t.params).filter(Boolean).flatMap(p => p.type))
 * ```
 * at runtime in the current codebase.
 */
export type Param = "text" | "number" | "multitext" | "text16" | "_multi_column:number" | "_multi_column:all";

//chatGPT to the rescue
export type BaseChartConstructor<TProps extends BaseConfig> = new (...args: any[]) => BaseChart<TProps>;
/**
 * Describes how a chart will be displayed in the 'add chart' dialog, with an optional `init` function that will
 * result in a sufficiently fully-formed configuration object for creating the chart.
 */
export type ChartType<T extends BaseConfig> = {
    /** A class extending BaseChart */
    class: BaseChartConstructor<T>;
    /** The human-readable name that will appear in the 'add chart' dialog etc. */
    name: string;

    /**
     * A list of columns that are required for the chart to function correctly, or a function that returns
     * a boolean based on evaluation of whether the data store meets arbitrary criteria defined by the chart.
     */
    required?: string[] | ((ds: DataStore) => boolean);
    /** 
     * If extra initialisation logic is required in order to change the internal configuration of the chart
     * into a form that can be used by the chart class, this function can be used to do that.
     * 
     * @param config A configuration object for the chart, in the intermediate form used internally by the Add Chart dialog.
     * @param dataSource The data source that the chart will be using (nb, do we mean `DataSource` or `DataStore`? - check, there are
     * frequent inconsistencies around this).
     */
    init?: (config: any, dataSource: any, extraControls: any) => void;
    /**
     * Specifies additional controls that will be displayed in the 'add chart' dialog, for configuring things other than what
     * goes into `params`.
     */
    extra_controls?: (dataStore: DataStore) => ExtraControl<GuiSpecType>[]; //not GuiSpec<GuiSpecType>[]; - AddChartDialog & SettingsDialog behave differently - would like to unify
    params?: { type: Param | Param[]; name: string }[];
    configEntriesUsingColumns?: string[];
    /** @deprecated in favour of {@link loadColumnData} decorator */
    methodsUsingColumns?: string[];
    allow_user_add?: boolean;
};

export type ChartTypeMap = Record<string, ChartType<any>>;

/** A dictionary of all the chart types, also accessible as `BaseChart.types` */
export const chartTypes: ChartTypeMap = {};
