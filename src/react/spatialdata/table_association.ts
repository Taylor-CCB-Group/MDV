import type { ShapesRenderData, SpatialData } from "@spatialdata/core";
import { loadAssociatedTableFeatureRows } from "@spatialdata/core";
import type { LayerConfig, LayerType, RenderStackLayerInputs } from "@spatialdata/vis";
import { useEffect, useMemo, useRef, useState } from "react";

import type DataStore from "@/datastore/DataStore";
import { useChartManager, useDataSources } from "@/react/hooks";

type ShapesLayerConfig = Extract<LayerConfig, { type: "shapes" }>;
type LabelsLayerConfig = Extract<LayerConfig, { type: "labels" }>;
type ShapeFeatureState = NonNullable<ShapesLayerConfig["featureState"]>;
type LabelFeatureState = NonNullable<LabelsLayerConfig["featureState"]>;
type AssociatedFeatureState = ShapeFeatureState | LabelFeatureState;
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
export type AssociableSpatialElementType = Extract<LayerType, "image" | "points" | "labels" | "shapes">;
type FillColorAssociableLayerType = Extract<AssociableSpatialElementType, "shapes" | "labels">;
type SpatialDataAssociationSource = {
    getAssociatedTables: (kind: SpatialDataAssociationKind, key: string) => Array<[string, unknown]>;
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

function spatialDataAssociationKind(elementType: AssociableSpatialElementType): SpatialDataAssociationKind {
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

function tableNamesForDataSource<TDataStore>(dataSource: DataSourceAssociationCandidate<TDataStore>): string[] {
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
    const rowIndexByFeatureId = new Map<string, number>();
    renderData.featureIds.forEach((featureId, featureIndex) => {
        const rowIndex = renderData.rowIndexByFeatureIndex[featureIndex];
        if (rowIndex === undefined) return;
        rowIndexByFeatureId.set(featureId, rowIndex);
    });
    return buildAssociatedFeatureStateFromRowMap({
        rowIndexByFeatureId,
        visibleRows,
        rowCount,
        baseFeatureState,
        colorForRow,
        alpha,
    });
}

export function buildAssociatedFeatureStateFromRowMap({
    rowIndexByFeatureId,
    visibleRows,
    rowCount,
    baseFeatureState,
    colorForRow,
    alpha,
}: {
    rowIndexByFeatureId: Map<string, number> | Iterable<[string, number]>;
    visibleRows: ArrayLike<number>;
    rowCount: number;
    baseFeatureState?: AssociatedFeatureState;
    colorForRow?: RowColorFunction;
    alpha: number;
}): AssociatedFeatureState | undefined {
    const visibleRowSet = new Set(Array.from(visibleRows));
    const hiddenFeatureIds = new Set(baseFeatureState?.hiddenFeatureIds ?? []);
    const fillColorByFeatureId = {
        ...(baseFeatureState?.fillColorByFeatureId ?? {}),
    };
    let hasTableState = false;

    for (const [featureId, rawRowIndex] of rowIndexByFeatureId) {
        const rowIndex = validAssociatedRow(rawRowIndex, rowCount);
        if (rowIndex === null) continue;

        hasTableState = true;
        if (!visibleRowSet.has(rowIndex)) {
            hiddenFeatureIds.add(featureId);
        }

        if (colorForRow) {
            const color = toRgba(colorForRow(rowIndex), alpha);
            if (color) fillColorByFeatureId[featureId] = color;
        }
    }

    if (!hasTableState && !baseFeatureState) return undefined;

    return {
        ...baseFeatureState,
        fillColorByFeatureId: Object.keys(fillColorByFeatureId).length > 0 ? fillColorByFeatureId : undefined,
        hiddenFeatureIds: hiddenFeatureIds.size > 0 ? Array.from(hiddenFeatureIds) : undefined,
    };
}

/**
 * Opt into upstream `fillColorByColumn` (skip MDV featureState fill colours).
 * Use with a linked `@spatialdata/vis` while verifying a fix before publish:
 *
 *   localStorage.MDV_USE_UPSTREAM_FILL_COLOR = "1"  // then reload
 */
export function preferUpstreamFillColorByColumn(): boolean {
    if (typeof window === "undefined") return false;
    try {
        return window.localStorage?.getItem("MDV_USE_UPSTREAM_FILL_COLOR") === "1";
    } catch {
        return false;
    }
}

/**
 * Upstream `@spatialdata/vis` still races `fillColorByColumn` loads against the canvas
 * for both shapes and labels (PR #119 last-good helps the load flash, but column switches
 * under MDV's adapter still fail to paint reliably). Drive colours through MDV
 * `featureState` instead and strip the column prop from viewer inputs so upstream does
 * not take that path — unless {@link preferUpstreamFillColorByColumn} is set.
 */
export function omitFillColorByColumn<T extends LayerConfig>(layer: T): T {
    if (!("fillColorByColumn" in layer)) return layer;
    const { fillColorByColumn: _removed, ...rest } = layer;
    return rest as T;
}

export function layerFillColorColumnName(layer: LayerConfig): string | undefined {
    if (!("fillColorByColumn" in layer)) return undefined;
    return layer.fillColorByColumn?.columnName;
}

/**
 * While a newly selected colour column is still loading, keep the previous feature colours
 * so the canvas does not flash through an unannotated/default state.
 */
export function withPreservedFillColorsWhileLoading({
    featureState,
    fillColumnName,
    colorReady,
    previousFillColorByFeatureId,
}: {
    featureState: AssociatedFeatureState | undefined;
    fillColumnName: string | undefined;
    colorReady: boolean;
    previousFillColorByFeatureId?: Record<string, RgbaColor>;
}): AssociatedFeatureState | undefined {
    if (!fillColumnName || colorReady || !previousFillColorByFeatureId) {
        return featureState;
    }
    return {
        ...featureState,
        fillColorByFeatureId: previousFillColorByFeatureId,
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

function labelElementKeys(layers: RenderStackLayerInputs["layers"], layerOrder: string[]) {
    const keys = new Set<string>();
    for (const layerId of layerOrder) {
        const layer = layers[layerId];
        if (layer?.type === "labels") keys.add(layer.elementKey);
    }
    return Array.from(keys);
}

export function useShapesRenderDataByElementKey(spatialData: SpatialData | undefined, elementKeys: string[]) {
    const [renderDataByElementKey, setRenderDataByElementKey] = useState<Record<string, ShapesRenderData | undefined>>(
        {},
    );

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
                    console.warn(`Failed to load SpatialData shapes render data for ${elementKey}`, error);
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
    }, [spatialData, elementKeys]);

    return renderDataByElementKey;
}

function useLabelsRowIndexByFeatureId(spatialData: SpatialData | undefined, elementKeys: string[]) {
    const [rowIndexByFeatureId, setRowIndexByFeatureId] = useState<Record<string, Map<string, number> | undefined>>({});

    useEffect(() => {
        let cancelled = false;
        if (!spatialData || elementKeys.length === 0) {
            setRowIndexByFeatureId({});
            return;
        }

        Promise.all(
            elementKeys.map(async (elementKey) => {
                try {
                    const rows = await loadAssociatedTableFeatureRows({
                        spatialData,
                        kind: "labels",
                        key: elementKey,
                    });
                    return [elementKey, rows.rowIndexByFeatureId] as const;
                } catch (error) {
                    console.warn(`Failed to load SpatialData labels table association for ${elementKey}`, error);
                    return [elementKey, undefined] as const;
                }
            }),
        ).then((entries) => {
            if (cancelled) return;
            setRowIndexByFeatureId(Object.fromEntries(entries));
        });

        return () => {
            cancelled = true;
        };
    }, [spatialData, elementKeys]);

    return rowIndexByFeatureId;
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
    return getShapesTableAssociation(renderData, table.dataStore.size, table.tableName, table.dataSourceName);
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
    }, [dataStores]);

    return version;
}

function isFillColorAssociableLayer(
    layer: LayerConfig | undefined,
): layer is Extract<LayerConfig, { type: FillColorAssociableLayerType }> {
    return layer?.type === "shapes" || layer?.type === "labels";
}

function getFillColumnsByDataSource(
    layers: RenderStackLayerInputs["layers"],
    layerOrder: string[],
    tableByAssociationKey: Record<string, AssociatedElementTable>,
) {
    const columnsByDataSource: Record<string, Set<string>> = {};
    for (const layerId of layerOrder) {
        const layer = layers[layerId];
        if (!isFillColorAssociableLayer(layer)) continue;
        const columnName = layerFillColorColumnName(layer);
        if (!columnName) continue;
        const table = tableByAssociationKey[`${layer.type}:${layer.elementKey}`];
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

function useLoadedColorColumnVersion(fillColumnsByDataSource: Record<string, string[]>) {
    const chartManager = useChartManager();
    const [version, setVersion] = useState(0);

    useEffect(() => {
        let cancelled = false;
        for (const [dataSourceName, columns] of Object.entries(fillColumnsByDataSource)) {
            const dataStore = chartManager.getDataSource(dataSourceName);
            const missing = columns.filter((columnName) => !dataStore.columnsWithData.includes(columnName));
            if (missing.length === 0) continue;
            chartManager.loadColumnSet(missing, dataSourceName, () => {
                if (!cancelled) setVersion((current) => current + 1);
            });
        }
        return () => {
            cancelled = true;
        };
    }, [chartManager, fillColumnsByDataSource]);

    return version;
}

function layerFillAlpha(layer: LayerConfig): number {
    if (layer.type === "shapes") return layer.fillColor?.[3] ?? 180;
    return 255;
}

export function useAssociatedShapesLayerInputs(
    spatialData: SpatialData | undefined,
    layerInputs: RenderStackLayerInputs,
): RenderStackLayerInputs {
    const dataSources = useDataSources();
    const shapeKeys = useMemo(
        () => shapeElementKeys(layerInputs.layers, layerInputs.layerOrder),
        [layerInputs.layers, layerInputs.layerOrder],
    );
    const labelKeys = useMemo(
        () => labelElementKeys(layerInputs.layers, layerInputs.layerOrder),
        [layerInputs.layers, layerInputs.layerOrder],
    );
    const renderDataByElementKey = useShapesRenderDataByElementKey(spatialData, shapeKeys);
    const labelsRowIndexByFeatureId = useLabelsRowIndexByFeatureId(spatialData, labelKeys);
    const tableByAssociationKey = useMemo(() => {
        const entries: Array<[string, AssociatedElementTable]> = [];
        for (const elementKey of shapeKeys) {
            entries.push([
                `shapes:${elementKey}`,
                resolveAssociatedElementTable({
                    spatialData,
                    elementType: "shapes",
                    elementKey,
                    dataSources,
                }),
            ]);
        }
        for (const elementKey of labelKeys) {
            entries.push([
                `labels:${elementKey}`,
                resolveAssociatedElementTable({
                    spatialData,
                    elementType: "labels",
                    elementKey,
                    dataSources,
                }),
            ]);
        }
        return Object.fromEntries(entries);
    }, [spatialData, shapeKeys, labelKeys, dataSources]);
    const fillColumnsByDataSource = useMemo(
        () => getFillColumnsByDataSource(layerInputs.layers, layerInputs.layerOrder, tableByAssociationKey),
        [layerInputs.layers, layerInputs.layerOrder, tableByAssociationKey],
    );
    const loadedColorColumnVersion = useLoadedColorColumnVersion(fillColumnsByDataSource);
    const associatedDataStores = useMemo(() => {
        const storesByName = new Map<string, DataStore>();
        for (const table of Object.values(tableByAssociationKey)) {
            if (table.status === "resolved") storesByName.set(table.dataSourceName, table.dataStore);
        }
        return Array.from(storesByName.values());
    }, [tableByAssociationKey]);
    const filterVersion = useDataStoreFilterVersion(associatedDataStores);
    const visibleRowsByDataSource = useMemo(() => {
        filterVersion;
        return Object.fromEntries(
            associatedDataStores.map((dataStore) => [dataStore.name, visibleRowsForDataStore(dataStore)]),
        );
    }, [associatedDataStores, filterVersion]);
    const colorFunctionByDataSource = useMemo(() => {
        loadedColorColumnVersion;
        return Object.fromEntries(
            Object.entries(fillColumnsByDataSource).map(([dataSourceName, columnNames]) => {
                const table = Object.values(tableByAssociationKey).find(
                    (candidate) => candidate.status === "resolved" && candidate.dataSourceName === dataSourceName,
                );
                const loadedColumnNames =
                    table?.status === "resolved" ? new Set(table.dataStore.columnsWithData) : new Set<string>();
                return [
                    dataSourceName,
                    table?.status === "resolved"
                        ? createColorFunctionByColumn(table.dataStore, columnNames, loadedColumnNames, {})
                        : {},
                ];
            }),
        );
    }, [fillColumnsByDataSource, loadedColorColumnVersion, tableByAssociationKey]);
    const lastFillColorsRef = useRef<Record<string, Record<string, RgbaColor>>>({});
    const useUpstreamFillColor = preferUpstreamFillColorByColumn();

    const layers = useMemo(() => {
        let changed = false;
        const nextLayers = { ...layerInputs.layers };

        for (const layerId of layerInputs.layerOrder) {
            const layer = layerInputs.layers[layerId];
            if (!isFillColorAssociableLayer(layer)) continue;

            const table = tableByAssociationKey[`${layer.type}:${layer.elementKey}`];
            if (table?.status !== "resolved") continue;

            const fillColumnName = layerFillColorColumnName(layer);
            // When verifying upstream fillColorByColumn, still project filter/hidden state
            // but leave fill colours to the viewer column path.
            const colorForRow =
                useUpstreamFillColor || !fillColumnName
                    ? undefined
                    : colorFunctionByDataSource[table.dataSourceName]?.[fillColumnName];
            const colorReady = useUpstreamFillColor || !fillColumnName || colorForRow !== undefined;
            const visibleRows =
                visibleRowsByDataSource[table.dataSourceName] ?? visibleRowsForDataStore(table.dataStore);

            let featureState: AssociatedFeatureState | undefined;
            if (layer.type === "shapes") {
                const renderData = renderDataByElementKey[layer.elementKey];
                if (!renderData) continue;
                featureState = buildAssociatedShapesFeatureState({
                    renderData,
                    visibleRows,
                    rowCount: table.dataStore.size,
                    baseFeatureState: layer.featureState,
                    colorForRow,
                    alpha: layerFillAlpha(layer),
                });
            } else {
                const rowIndexByFeatureId = labelsRowIndexByFeatureId[layer.elementKey];
                if (!rowIndexByFeatureId) continue;
                featureState = buildAssociatedFeatureStateFromRowMap({
                    rowIndexByFeatureId,
                    visibleRows,
                    rowCount: table.dataStore.size,
                    baseFeatureState: layer.featureState,
                    colorForRow,
                    alpha: layerFillAlpha(layer),
                });
            }

            if (!useUpstreamFillColor) {
                featureState = withPreservedFillColorsWhileLoading({
                    featureState,
                    fillColumnName,
                    colorReady,
                    previousFillColorByFeatureId: lastFillColorsRef.current[layerId],
                });

                if (!fillColumnName) {
                    const { [layerId]: _removed, ...remaining } = lastFillColorsRef.current;
                    lastFillColorsRef.current = remaining;
                } else if (colorReady && featureState?.fillColorByFeatureId) {
                    lastFillColorsRef.current[layerId] = featureState.fillColorByFeatureId;
                }
            }

            if (!featureState && !fillColumnName) continue;

            const projected = useUpstreamFillColor ? layer : omitFillColorByColumn(layer);
            nextLayers[layerId] = featureState ? { ...projected, featureState } : projected;
            changed = true;
        }

        return changed ? nextLayers : layerInputs.layers;
    }, [
        layerInputs.layers,
        layerInputs.layerOrder,
        renderDataByElementKey,
        labelsRowIndexByFeatureId,
        tableByAssociationKey,
        visibleRowsByDataSource,
        colorFunctionByDataSource,
        useUpstreamFillColor,
    ]);

    return useMemo(
        () => ({
            layers,
            layerOrder: layerInputs.layerOrder,
        }),
        [layers, layerInputs.layerOrder],
    );
}
