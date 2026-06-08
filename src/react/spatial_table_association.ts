import type { SpatialData } from "@spatialdata/core";
import type DataStore from "@/datastore/DataStore";

export type TableAssociationStatus = "resolved" | "ambiguous" | "none";

export type TableAssociation = {
    status: TableAssociationStatus;
    tableName?: string;
    elementKey?: string;
    instanceKeyColumn?: string;
    spatialdataPath?: string;
    candidates?: string[];
};

type InferTableAssociationArgs = {
    spatialData: SpatialData | null | undefined;
    spatialdataPath: string | null | undefined;
    regionId: string | null | undefined;
    elementKey?: string | null;
    dataStore: DataStore;
};

function parseElementKeyFromRegion(regionId: string, spatialdataPath: string): string | null {
    const prefix = `${spatialdataPath}_`;
    if (!regionId.startsWith(prefix)) return null;
    return regionId.slice(prefix.length) || null;
}

function tableNamesFromMdvRows(
    dataStore: DataStore,
    spatialdataPath: string,
    regionId: string,
): string[] {
    const tableNameColumn = dataStore.columnIndex.table_name;
    const spatialPathColumn = dataStore.columnIndex.spatialdata_path;
    const regionColumn = dataStore.columnIndex.spatial_region;
    if (!tableNameColumn || !spatialPathColumn || !regionColumn) {
        return [];
    }

    const names = new Set<string>();
    const pathValues = spatialPathColumn.values;
    const regionValues = regionColumn.values;
    const tableValues = tableNameColumn.values;
    if (!pathValues || !regionValues || !tableValues) return [];

    const length = Math.min(pathValues.length, regionValues.length, tableValues.length);
    for (let i = 0; i < length; i += 1) {
        if (pathValues[i] === spatialdataPath && regionValues[i] === regionId) {
            const name = tableValues[i];
            if (typeof name === "string" && name.length > 0) {
                names.add(name);
            }
        }
    }
    return [...names];
}

function tableNamesFromSpatialData(
    spatialData: SpatialData,
    elementKey: string,
): string[] {
    const tables = spatialData.tables ?? {};
    const matches: string[] = [];

    for (const [tableName, table] of Object.entries(tables)) {
        const tableRecord = table as { uns?: { spatialdata_attrs?: { region?: string | string[] } } };
        const attrs = tableRecord.uns?.spatialdata_attrs;
        const region = attrs?.region;
        if (!region) continue;
        const regions = Array.isArray(region) ? region : [region];
        if (regions.includes(elementKey)) {
            matches.push(tableName);
        }
    }

    return matches;
}

function instanceKeyColumnForTable(
    spatialData: SpatialData,
    tableName: string,
    dataStore: DataStore,
): string | undefined {
    const table = spatialData.tables?.[tableName] as
        | { uns?: { spatialdata_attrs?: { instance_key?: string } } }
        | undefined;
    const attrs = table?.uns?.spatialdata_attrs;
    const instanceKey = attrs?.instance_key;
    if (typeof instanceKey === "string" && dataStore.columnIndex[instanceKey]) {
        return instanceKey;
    }
    return undefined;
}

function reconcileTableNames(
    mdvTables: string[],
    spatialTables: string[],
    allSpatialTableNames: string[],
): { status: TableAssociationStatus; tableName?: string; candidates?: string[] } {
    if (allSpatialTableNames.length === 1) {
        return { status: "resolved", tableName: allSpatialTableNames[0] };
    }

    const mdvSet = new Set(mdvTables);
    const spatialSet = new Set(spatialTables);
    const intersection = [...mdvSet].filter((name) => spatialSet.has(name));

    if (intersection.length === 1) {
        return { status: "resolved", tableName: intersection[0] };
    }

    const union = [...new Set([...mdvTables, ...spatialTables])];
    if (union.length === 1) {
        return { status: "resolved", tableName: union[0] };
    }

    if (union.length > 1) {
        return { status: "ambiguous", candidates: union };
    }

    return { status: "none" };
}

export function inferTableAssociation({
    spatialData,
    spatialdataPath,
    regionId,
    elementKey,
    dataStore,
}: InferTableAssociationArgs): TableAssociation {
    if (!spatialdataPath) {
        return { status: "none" };
    }

    const resolvedElementKey =
        elementKey ??
        (regionId ? parseElementKeyFromRegion(regionId, spatialdataPath) : null);

    const allSpatialTableNames = spatialData ? Object.keys(spatialData.tables ?? {}) : [];
    if (allSpatialTableNames.length === 1) {
        const tableName = allSpatialTableNames[0];
        return {
            status: "resolved",
            tableName,
            elementKey: resolvedElementKey ?? undefined,
            instanceKeyColumn:
                spatialData ? instanceKeyColumnForTable(spatialData, tableName, dataStore) : undefined,
            spatialdataPath,
        };
    }

    const mdvTables =
        regionId ? tableNamesFromMdvRows(dataStore, spatialdataPath, regionId) : [];
    const spatialTables =
        spatialData && resolvedElementKey
            ? tableNamesFromSpatialData(spatialData, resolvedElementKey)
            : [];

    const reconciliation = reconcileTableNames(
        mdvTables,
        spatialTables,
        allSpatialTableNames,
    );

    const tableName = reconciliation.tableName;
    return {
        status: reconciliation.status,
        tableName,
        elementKey: resolvedElementKey ?? undefined,
        instanceKeyColumn:
            spatialData && tableName
                ? instanceKeyColumnForTable(spatialData, tableName, dataStore)
                : undefined,
        spatialdataPath,
        candidates: reconciliation.candidates,
    };
}
