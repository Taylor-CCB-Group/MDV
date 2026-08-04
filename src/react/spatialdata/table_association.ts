import type { SpatialData, ShapesRenderData } from "@spatialdata/core";
import type { LayerConfig, LayerType, RenderStackLayerInputs } from "@spatialdata/vis";
import { useEffect, useMemo, useState } from "react";

import type DataStore from "@/datastore/DataStore";
import { useChartManager, useDataSources } from "@/react/hooks";

type ShapesLayerConfig = Extract<LayerConfig, { type: "shapes" }>;
type ShapeFeatureState = NonNullable<ShapesLayerConfig["featureState"]>;
type RgbColor = [number, number, number];
type RgbaColor = [number, number, number, number];
type RowColorFunction = (rowIndex: number) => RgbColor | RgbaColor | undefined;
type ShapesRenderDataEntry = [string, ShapesRenderData | undefined];
type SpatialDataTableMetadata = { table_name?: unknown; table_id?: unknown };
type SpatialDataTablesMetadata = { tables?: SpatialDataTableMetadata[] };
type DataSourceColumnMetadata = {
    field?: unknown;
    name?: unknown;
    values?: unknown;
};
type DataSourceAssociationConfig = {
    columns?: DataSourceColumnMetadata[];
    spatialdata_tables?: SpatialDataTablesMetadata;
};
type SpatialDataAssociationKind = "images" | "points" | "labels" | "shapes";
export type AssociableSpatialElementType = Extract<
    LayerType,
    "image" | "points" | "labels" | "shapes"
>;
type SpatialDataAssociationSource = {
    getAssociatedTables: (
        kind: SpatialDataAssociationKind,
        key: string,
    ) => Array<[string, unknown]>;
};
type SpatialDataAwareDataStore = DataStore & {
    config: DataStore["config"] & { spatialdata_tables?: SpatialDataTablesMetadata };
};
export type DataSourceAssociationCandidate<TDataStore = Record<string, unknown>> = {
    name: string;
    dataStore: TDataStore & { config?: DataSourceAssociationConfig };
};
export type AssociatedDataSource = DataSourceAssociationCandidate<SpatialDataAwareDataStore>;
type AssociatedElementTable<TDataStore = SpatialDataAwareDataStore> =
    | { status: "none" }
    | { status: "ambiguous"; tableNames: string[]; dataSourceNames?: string[] }
    | {
          status: "resolved";
          tableName: string;
          dataSourceName: string;
          dataStore: TDataStore;
      };

export type TableAssociation =
    | { status: "none" }
    | { status: "loading" }
    | { status: "ambiguous" }
    | {
          status: "resolved";
          tableName: string;
          dataSourceName: string;
          matchedFeatureCount?: number;
          featureCount?: number;
      };

export const NO_TABLE_ASSOCIATION: TableAssociation = { status: "none" };

function validAssociatedRow(rowIndex: number | undefined, rowCount: number): number | null {
    if (rowIndex === undefined || rowIndex < 0 || rowIndex >= rowCount) return null;
    return rowIndex;
}

export function getShapesTableAssociation(
    renderData: ShapesRenderData | undefined,
    rowCount: number,
    tableName: string,
    dataSourceName: string,
): TableAssociation {
    if (!renderData) return NO_TABLE_ASSOCIATION;

    let matchedFeatureCount = 0;
    for (const rowIndex of renderData.rowIndexByFeatureIndex) {
        if (validAssociatedRow(rowIndex, rowCount) !== null) matchedFeatureCount++;
    }

    if (matchedFeatureCount === 0) return NO_TABLE_ASSOCIATION;
    return {
        status: "resolved",
        tableName,
        dataSourceName,
        matchedFeatureCount,
        featureCount: renderData.featureIds.length,
    };
}

function associatedSpatialDataTableNames(
    spatialData: SpatialDataAssociationSource | undefined,
    elementType: AssociableSpatialElementType,
    elementKey: string | undefined,
): string[] {
    if (!spatialData || !elementKey) return [];
    return spatialData
        .getAssociatedTables(spatialDataAssociationKind(elementType), elementKey)
        .map(([tableName]) => tableName);
}

function spatialDataAssociationKind(
    elementType: AssociableSpatialElementType,
): SpatialDataAssociationKind {
    switch (elementType) {
        case "image":
            return "images";
        case "labels":
            return "labels";
        case "points":
            return "points";
        case "shapes":
            return "shapes";
    }
}

function tableNamesForDataSource<TDataStore>(
    dataSource: DataSourceAssociationCandidate<TDataStore>,
): string[] {
    const names = new Set<string>();
    const provenance = dataSource.dataStore.config?.spatialdata_tables;
    const tables = Array.isArray(provenance?.tables) ? provenance.tables : [];
    for (const table of tables) {
        if (typeof table.table_name === "string") names.add(table.table_name);
        if (typeof table.table_id === "string") {
            const tableIdLeaf = table.table_id.split("/").pop();
            if (tableIdLeaf) names.add(tableIdLeaf);
        }
    }

    const columns = dataSource.dataStore.config?.columns;
    if (Array.isArray(columns)) {
        for (const column of columns) {
            const columnId =
                typeof column.field === "string"
                    ? column.field
                    : typeof column.name === "string"
                      ? column.name
                      : undefined;
            if (columnId !== "table_name") continue;
            if (!Array.isArray(column.values)) continue;
            for (const value of column.values) {
                if (typeof value === "string") names.add(value);
            }
        }
    }

    return Array.from(names);
}

export function resolveAssociatedElementTable<TDataStore>({
    spatialData,
    elementType,
    elementKey,
    dataSources,
}: {
    spatialData: SpatialDataAssociationSource | undefined;
    elementType: AssociableSpatialElementType;
    elementKey: string | undefined;
    dataSources: DataSourceAssociationCandidate<TDataStore>[];
}): AssociatedElementTable<TDataStore> {
    const tableNames = associatedSpatialDataTableNames(spatialData, elementType, elementKey);
    if (tableNames.length === 0) return { status: "none" };
    if (tableNames.length > 1) return { status: "ambiguous", tableNames };

    const tableName = tableNames[0];
    const matchingDataSources = dataSources.filter((dataSource) => {
        if (dataSource.name === tableName) return true;
        return tableNamesForDataSource(dataSource).includes(tableName);
    });

    if (matchingDataSources.length === 0) return { status: "none" };
    if (matchingDataSources.length > 1) {
        return {
            status: "ambiguous",
            tableNames,
            dataSourceNames: matchingDataSources.map((dataSource) => dataSource.name),
        };
    }

    const dataSource = matchingDataSources[0];
    return {
        status: "resolved",
        tableName,
        dataSourceName: dataSource.name,
        dataStore: dataSource.dataStore,
    };
}

function toRgba(color: RgbColor | RgbaColor | undefined, alpha: number): RgbaColor | null {
    if (!color) return null;
    return [
        Math.max(0, Math.min(255, color[0] ?? 0)),
        Math.max(0, Math.min(255, color[1] ?? 0)),
        Math.max(0, Math.min(255, color[2] ?? 0)),
        Math.max(0, Math.min(255, color[3] ?? alpha)),
    ];
}

export function buildAssociatedShapesFeatureState({
    renderData,
    visibleRows,
    rowCount,
    baseFeatureState,
    colorForRow,
    alpha,
}: {
    renderData: ShapesRenderData;
    visibleRows: ArrayLike<number>;
    rowCount: number;
    baseFeatureState?: ShapeFeatureState;
    colorForRow?: RowColorFunction;
    alpha: number;
}): ShapeFeatureState | undefined {
    const visibleRowSet = new Set(Array.from(visibleRows));
    const hiddenFeatureIds = new Set(baseFeatureState?.hiddenFeatureIds ?? []);
    const fillColorByFeatureId = {
        ...(baseFeatureState?.fillColorByFeatureId ?? {}),
    };
    let hasTableState = false;

    renderData.featureIds.forEach((featureId, featureIndex) => {
        const rowIndex = validAssociatedRow(
            renderData.rowIndexByFeatureIndex[featureIndex],
            rowCount,
        );
        if (rowIndex === null) return;

        hasTableState = true;
        if (!visibleRowSet.has(rowIndex)) {
            hiddenFeatureIds.add(featureId);
        }

        if (colorForRow) {
            const color = toRgba(colorForRow(rowIndex), alpha);
            if (color) fillColorByFeatureId[featureId] = color;
        }
    });

    if (!hasTableState && !baseFeatureState) return undefined;

    return {
        ...baseFeatureState,
        fillColorByFeatureId:
            Object.keys(fillColorByFeatureId).length > 0
                ? fillColorByFeatureId
                : undefined,
        hiddenFeatureIds:
            hiddenFeatureIds.size > 0 ? Array.from(hiddenFeatureIds) : undefined,
    };
}

function shapeElementKeys(layers: RenderStackLayerInputs["layers"], layerOrder: string[]) {
    const keys = new Set<string>();
    for (const layerId of layerOrder) {
        const layer = layers[layerId];
        if (layer?.type === "shapes") keys.add(layer.elementKey);
    }
    return Array.from(keys);
}

export function useShapesRenderDataByElementKey(
    spatialData: SpatialData | undefined,
    elementKeys: string[],
) {
    const keyFingerprint = elementKeys.join("\0");
    const [renderDataByElementKey, setRenderDataByElementKey] = useState<
        Record<string, ShapesRenderData | undefined>
    >({});

    useEffect(() => {
        let cancelled = false;
        if (!spatialData || elementKeys.length === 0) {
            setRenderDataByElementKey({});
            return;
        }

        Promise.all(
            elementKeys.map(async (elementKey): Promise<ShapesRenderDataEntry> => {
                const element = spatialData.shapes?.[elementKey];
                if (!element) return [elementKey, undefined];
                try {
                    return [elementKey, await element.loadRenderData()];
                } catch (error) {
                    console.warn(
                        `Failed to load SpatialData shapes render data for ${elementKey}`,
                        error,
                    );
                    return [elementKey, undefined];
                }
            }),
        ).then((entries) => {
            if (cancelled) return;
            setRenderDataByElementKey(Object.fromEntries(entries));
        });

        return () => {
            cancelled = true;
        };
    }, [spatialData, keyFingerprint]);

    return renderDataByElementKey;
}

export function useElementTableAssociation(
    spatialData: SpatialData | undefined,
    elementType: AssociableSpatialElementType,
    elementKey: string | undefined,
    dataSources: AssociatedDataSource[],
): TableAssociation {
    const elementKeys = useMemo(
        () => (elementType === "shapes" && elementKey ? [elementKey] : []),
        [elementType, elementKey],
    );
    const renderDataByElementKey = useShapesRenderDataByElementKey(spatialData, elementKeys);
    if (!elementKey || !spatialData) return NO_TABLE_ASSOCIATION;
    const table = resolveAssociatedElementTable({
        spatialData,
        elementType,
        elementKey,
        dataSources,
    });
    if (table.status === "ambiguous") return { status: "ambiguous" };
    if (table.status !== "resolved") return NO_TABLE_ASSOCIATION;
    if (elementType !== "shapes") {
        return {
            status: "resolved",
            tableName: table.tableName,
            dataSourceName: table.dataSourceName,
        };
    }
    const renderData = renderDataByElementKey[elementKey];
    if (!renderData) return { status: "loading" };
    return getShapesTableAssociation(
        renderData,
        table.dataStore.size,
        table.tableName,
        table.dataSourceName,
    );
}

export function useShapesTableAssociation(
    spatialData: SpatialData | undefined,
    elementKey: string | undefined,
    dataSources: AssociatedDataSource[],
): TableAssociation {
    return useElementTableAssociation(spatialData, "shapes", elementKey, dataSources);
}

function createColorFunctionByColumn(
    dataStore: DataStore,
    columnNames: string[],
    loadedColumnNames: Set<string>,
    colorOptions: {
        log_color_scale?: boolean;
        fallbackOnZero?: boolean;
        hideMissing?: boolean;
    },
) {
    const colorFunctionByColumn: Record<string, RowColorFunction | undefined> = {};
    for (const columnName of columnNames) {
        if (!loadedColumnNames.has(columnName)) continue;
        colorFunctionByColumn[columnName] = dataStore.getColorFunction(columnName, {
            asArray: true,
            overideValues: {
                colorLogScale: colorOptions.log_color_scale,
                fallbackOnZero: colorOptions.fallbackOnZero,
                hideMissing: colorOptions.hideMissing,
            },
        });
    }
    return colorFunctionByColumn;
}

function visibleRowsForDataStore(dataStore: DataStore): Uint32Array {
    const visibleRows = new Uint32Array(dataStore.filterSize);
    let index = 0;
    for (let rowIndex = 0; rowIndex < dataStore.size; rowIndex++) {
        if (dataStore.filterArray[rowIndex] !== 0) continue;
        visibleRows[index++] = rowIndex;
    }
    return visibleRows;
}

function useDataStoreFilterVersion(dataStores: DataStore[]) {
    const dataStoreNames = dataStores.map((dataStore) => dataStore.name).join("\0");
    const [version, setVersion] = useState(0);

    useEffect(() => {
        const listenerId = `spatial-table-association-${Math.random().toString(36).slice(2)}`;
        const listener = (type: string) => {
            if (type === "filtered" || type === "data_added") {
                setVersion((current) => current + 1);
            }
        };

        for (const dataStore of dataStores) {
            dataStore.addListener(listenerId, listener);
        }
        return () => {
            for (const dataStore of dataStores) {
                dataStore.removeListener(listenerId);
            }
        };
    }, [dataStores, dataStoreNames]);

    return version;
}

function getFillColumnsByDataSource(
    layers: RenderStackLayerInputs["layers"],
    layerOrder: string[],
    tableByElementKey: Record<string, AssociatedElementTable>,
) {
    const columnsByDataSource: Record<string, Set<string>> = {};
    for (const layerId of layerOrder) {
        const layer = layers[layerId];
        if (layer?.type !== "shapes") continue;
        const columnName = layer.fillColorByColumn?.columnName;
        if (!columnName) continue;
        const table = tableByElementKey[layer.elementKey];
        if (table?.status !== "resolved") continue;
        columnsByDataSource[table.dataSourceName] ??= new Set();
        columnsByDataSource[table.dataSourceName]?.add(columnName);
    }
    return Object.fromEntries(
        Object.entries(columnsByDataSource).map(([dataSourceName, columns]) => [
            dataSourceName,
            Array.from(columns).sort(),
        ]),
    );
}

function useLoadedColorColumnVersion(
    fillColumnsByDataSource: Record<string, string[]>,
) {
    const chartManager = useChartManager();
    const columnFingerprint = Object.entries(fillColumnsByDataSource)
        .map(([dataSourceName, columns]) => `${dataSourceName}:${columns.join(",")}`)
        .join("|");
    const [version, setVersion] = useState(0);

    useEffect(() => {
        let cancelled = false;
        for (const [dataSourceName, columns] of Object.entries(fillColumnsByDataSource)) {
            const dataStore = chartManager.getDataSource(dataSourceName);
            const missing = columns.filter(
                (columnName) => !dataStore.columnsWithData.includes(columnName),
            );
            if (missing.length === 0) continue;
            chartManager.loadColumnSet(missing, dataSourceName, () => {
                if (!cancelled) setVersion((current) => current + 1);
            });
        }
        return () => {
            cancelled = true;
        };
    }, [chartManager, columnFingerprint, fillColumnsByDataSource]);

    return version;
}

export function useAssociatedShapesLayerInputs(
    spatialData: SpatialData | undefined,
    layerInputs: RenderStackLayerInputs,
): RenderStackLayerInputs {
    const dataSources = useDataSources();
    const elementKeys = useMemo(
        () => shapeElementKeys(layerInputs.layers, layerInputs.layerOrder),
        [layerInputs.layers, layerInputs.layerOrder],
    );
    const renderDataByElementKey = useShapesRenderDataByElementKey(spatialData, elementKeys);
    const tableByElementKey = useMemo(
        () =>
            Object.fromEntries(
                elementKeys.map((elementKey) => [
                    elementKey,
                    resolveAssociatedElementTable({
                        spatialData,
                        elementType: "shapes",
                        elementKey,
                        dataSources,
                    }),
                ]),
            ),
        [spatialData, elementKeys, dataSources],
    );
    const fillColumnsByDataSource = useMemo(
        () =>
            getFillColumnsByDataSource(
                layerInputs.layers,
                layerInputs.layerOrder,
                tableByElementKey,
            ),
        [layerInputs.layers, layerInputs.layerOrder, tableByElementKey],
    );
    const loadedColorColumnVersion =
        useLoadedColorColumnVersion(fillColumnsByDataSource);
    const associatedDataStores = useMemo(() => {
        const storesByName = new Map<string, DataStore>();
        for (const table of Object.values(tableByElementKey)) {
            if (table.status === "resolved") storesByName.set(table.dataSourceName, table.dataStore);
        }
        return Array.from(storesByName.values());
    }, [tableByElementKey]);
    const filterVersion = useDataStoreFilterVersion(associatedDataStores);
    const visibleRowsByDataSource = useMemo(() => {
        filterVersion;
        return Object.fromEntries(
            associatedDataStores.map((dataStore) => [
                dataStore.name,
                visibleRowsForDataStore(dataStore),
            ]),
        );
    }, [associatedDataStores, filterVersion]);
    const colorFunctionByDataSource = useMemo(
        () => {
            loadedColorColumnVersion;
            return Object.fromEntries(
                Object.entries(fillColumnsByDataSource).map(([dataSourceName, columnNames]) => {
                    const table = Object.values(tableByElementKey).find(
                        (candidate) =>
                            candidate.status === "resolved" &&
                            candidate.dataSourceName === dataSourceName,
                    );
                    const loadedColumnNames =
                        table?.status === "resolved"
                            ? new Set(table.dataStore.columnsWithData)
                            : new Set<string>();
                    return [
                        dataSourceName,
                        table?.status === "resolved"
                            ? createColorFunctionByColumn(
                                  table.dataStore,
                                  columnNames,
                                  loadedColumnNames,
                                  {},
                              )
                            : {},
                    ];
                }),
            );
        },
        [
            fillColumnsByDataSource,
            loadedColorColumnVersion,
            tableByElementKey,
        ],
    );

    const layers = useMemo(
        () => {
            let changed = false;
            const nextLayers = { ...layerInputs.layers };

            for (const layerId of layerInputs.layerOrder) {
                const layer = layerInputs.layers[layerId];
                if (layer?.type !== "shapes") continue;

                const table = tableByElementKey[layer.elementKey];
                if (table?.status !== "resolved") continue;
                const renderData = renderDataByElementKey[layer.elementKey];
                if (!renderData) continue;

                const fillColumnName = layer.fillColorByColumn?.columnName;
                const featureState = buildAssociatedShapesFeatureState({
                    renderData,
                    visibleRows:
                        visibleRowsByDataSource[table.dataSourceName] ??
                        visibleRowsForDataStore(table.dataStore),
                    rowCount: table.dataStore.size,
                    baseFeatureState: layer.featureState,
                    colorForRow: fillColumnName
                        ? colorFunctionByDataSource[table.dataSourceName]?.[fillColumnName]
                        : undefined,
                    alpha: layer.fillColor?.[3] ?? 180,
                });

                if (!featureState) continue;
                nextLayers[layerId] = { ...layer, featureState };
                changed = true;
            }

            return changed ? nextLayers : layerInputs.layers;
        },
        [
            layerInputs.layers,
            layerInputs.layerOrder,
            renderDataByElementKey,
            tableByElementKey,
            visibleRowsByDataSource,
            colorFunctionByDataSource,
        ],
    );

    return useMemo(
        () => ({
            layers,
            layerOrder: layerInputs.layerOrder,
        }),
        [layers, layerInputs.layerOrder],
    );
}
