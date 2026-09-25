/**
 * `(region, instance id)` → DataStore row, for a SpatialData table converted into an MDV
 * datasource.
 *
 * `@spatialdata/core` aligns features to rows of ONE zarr table. The converter may have
 * concatenated several tables into one datasource (`table_handling` `by-region` or
 * `merge`), so those rows are only DataStore rows when the datasource is exactly that
 * table. Otherwise we join on values the converter already wrote to every row:
 * `spatialdata_table_id`, plus the table's own region-key and instance-key obs columns.
 *
 * The join key is always `(region, instance id)`: instance ids are unique only within a
 * region, which is what `region_key` exists to say.
 */

/** The part of a DataStore this needs; a real `DataStore` satisfies it. */
export type TableRowSource = {
    size: number;
    columnsWithData: string[];
    columnIndex: Record<string, { getValue: (index: number) => string | number } | undefined>;
};

export type TableRowPlan =
    /** The datasource is the zarr table, in its own order. */
    | { kind: "positional" }
    | {
          kind: "value";
          tableId: string;
          /** Fixed region for every row of the table, when it annotates only one. */
          region?: string;
          /** Column holding each row's region, when the table annotates several. */
          regionColumn?: string;
          instanceColumn: string;
          /** Every column the join reads; they must be loaded first. */
          columns: string[];
      }
    | { kind: "unresolved"; reason: string };

export type TableRowResolver = {
    /**
     * DataStore row for each instance id in `region`; `-1` where there is none.
     *
     * `tableRows` are the zarr-table rows `@spatialdata/core` aligned to the same ids.
     * Only the positional path reads them — the value join ignores them.
     */
    resolveRows(region: string, ids: ArrayLike<string>, tableRows: ArrayLike<number>): Int32Array;
};

type TableProvenance = {
    spatialdata_name?: unknown;
    table_name?: unknown;
    table_id?: unknown;
    region?: unknown;
    region_key?: unknown;
    instance_key?: unknown;
};

const TABLE_ID_COLUMN = "spatialdata_table_id";

function provenanceTables(provenance: unknown): TableProvenance[] {
    if (!provenance || typeof provenance !== "object" || !("tables" in provenance)) return [];
    return Array.isArray(provenance.tables) ? provenance.tables : [];
}

function regionsOf(table: TableProvenance): string[] {
    if (typeof table.region === "string") return [table.region];
    if (Array.isArray(table.region)) return table.region.filter((region) => typeof region === "string");
    return [];
}

/**
 * Decide how rows of `tableName` (from the SpatialData object `spatialDataName`) reach
 * rows of a datasource, from the `spatialdata_tables` provenance the converter wrote
 * into its config.
 */
export function planTableRows(
    provenance: unknown,
    tableName: string,
    spatialDataName: string | undefined,
): TableRowPlan {
    const tables = provenanceTables(provenance);
    // A single-table group keeps the table's own row order. No provenance at all is a
    // datasource matched by name, from before the converter recorded any.
    if (tables.length <= 1) return { kind: "positional" };

    const named = tables.filter((table) => table.table_name === tableName);
    const matches =
        spatialDataName === undefined ? named : named.filter((table) => table.spatialdata_name === spatialDataName);
    if (matches.length !== 1) {
        return {
            kind: "unresolved",
            reason: `${matches.length} tables named '${tableName}' from '${spatialDataName ?? "?"}'`,
        };
    }

    const table = matches[0];
    if (typeof table.table_id !== "string" || typeof table.instance_key !== "string" || !table.instance_key) {
        return { kind: "unresolved", reason: `table '${tableName}' has no table_id or instance_key` };
    }
    const regions = regionsOf(table);
    const regionColumn = regions.length > 1 && typeof table.region_key === "string" ? table.region_key : undefined;
    if (regions.length !== 1 && !regionColumn) {
        return { kind: "unresolved", reason: `table '${tableName}' has no usable region` };
    }
    return {
        kind: "value",
        tableId: table.table_id,
        region: regionColumn ? undefined : regions[0],
        regionColumn,
        instanceColumn: table.instance_key,
        columns: [TABLE_ID_COLUMN, ...(regionColumn ? [regionColumn] : []), table.instance_key],
    };
}

function positionalResolver(size: number): TableRowResolver {
    return {
        resolveRows(_region, ids, tableRows) {
            const rows = new Int32Array(ids.length).fill(-1);
            for (let i = 0; i < ids.length; i++) {
                const row = tableRows[i];
                if (row !== undefined && row >= 0 && row < size) rows[i] = row;
            }
            return rows;
        },
    };
}

/** `\0` cannot occur in a zarr element name. */
function joinKey(region: string, id: string) {
    return `${region}\0${id}`;
}

function valueResolver(source: TableRowSource, plan: Extract<TableRowPlan, { kind: "value" }>): TableRowResolver {
    const tableIds = source.columnIndex[TABLE_ID_COLUMN];
    const instances = source.columnIndex[plan.instanceColumn];
    const regions = plan.regionColumn ? source.columnIndex[plan.regionColumn] : undefined;
    const rowByKey = new Map<string, number>();
    if (tableIds && instances) {
        for (let row = 0; row < source.size; row++) {
            if (tableIds.getValue(row) !== plan.tableId) continue;
            const region = regions ? String(regions.getValue(row)) : plan.region;
            if (region === undefined) continue;
            const key = joinKey(region, String(instances.getValue(row)));
            if (!rowByKey.has(key)) rowByKey.set(key, row);
        }
    }
    return {
        resolveRows(region, ids) {
            const rows = new Int32Array(ids.length);
            for (let i = 0; i < ids.length; i++) {
                rows[i] = rowByKey.get(joinKey(region, ids[i])) ?? -1;
            }
            return rows;
        },
    };
}

/**
 * Built once per datasource, size and plan, and reused by every element and layer that
 * asks — including a batch of unique ids from a future points → cells join. Keyed by
 * value so a fresh-but-equal plan does not rebuild the join index.
 */
const resolverCache = new WeakMap<TableRowSource, Map<string, TableRowResolver>>();

/**
 * The resolver for `plan` over `source`, or `undefined` while the columns the join reads
 * have not loaded (or the plan cannot be resolved at all).
 */
export function getTableRowResolver(source: TableRowSource, plan: TableRowPlan): TableRowResolver | undefined {
    if (plan.kind === "unresolved") return undefined;
    if (plan.kind === "value" && !plan.columns.every((column) => source.columnsWithData.includes(column))) {
        return undefined;
    }

    const key =
        plan.kind === "positional"
            ? `positional\0${source.size}`
            : [source.size, plan.tableId, plan.region ?? "", plan.regionColumn ?? "", plan.instanceColumn].join("\0");
    let bySource = resolverCache.get(source);
    if (!bySource) {
        bySource = new Map();
        resolverCache.set(source, bySource);
    }
    let resolver = bySource.get(key);
    if (!resolver) {
        resolver = plan.kind === "positional" ? positionalResolver(source.size) : valueResolver(source, plan);
        bySource.set(key, resolver);
    }
    return resolver;
}

/** Feature id → DataStore row, dropping features with no row. */
export function datasourceRowsByFeatureId(
    resolver: TableRowResolver,
    region: string,
    featureIds: ArrayLike<string>,
    tableRows: ArrayLike<number>,
): Map<string, number> {
    const rows = resolver.resolveRows(region, featureIds, tableRows);
    const rowByFeatureId = new Map<string, number>();
    for (let i = 0; i < rows.length; i++) {
        if (rows[i] >= 0) rowByFeatureId.set(featureIds[i], rows[i]);
    }
    return rowByFeatureId;
}
