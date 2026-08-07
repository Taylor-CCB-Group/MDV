import type { ShapesRenderData, SpatialData } from "@spatialdata/core";
import { loadAssociatedTableFeatureRows } from "@spatialdata/core";
import type { LayerConfig, LayerType, RenderStackLayerInputs } from "@spatialdata/vis";
import { useEffect, useMemo, useRef, useState } from "react";

import type DataStore from "@/datastore/DataStore";
import { useChartManager, useDataSources } from "@/react/hooks";
import { fillColorSchemeFromDataStore } from "@/react/spatialdata/fill_color_scheme";

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
type AssociatedTableElement = { getObsColumnNames?: () => string[] };
type SpatialDataAssociationSource = {
    getAssociatedTables: (kind: SpatialDataAssociationKind, key: string) => Array<[string, AssociatedTableElement]>;
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
 * The obs columns of the table annotating an element — the exact set the viewer is
 * able to colour by on its own.
 *
 * This is the line the fill-colour routing is drawn on. MDV's column picker offers
 * everything in the DataStore, which is a SUPERSET of the table's obs: linked gene
 * scores, `mdv_cell_id`, anything computed at runtime. The viewer cannot read those
 * at all, so they stay on the per-feature path; everything in obs is handed over as
 * `fillColorByColumn` instead, with MDV's palette attached.
 *
 * Synchronous: opening the store already reads every node's attributes, so this is
 * a question we can answer before deciding how to colour, not after a load.
 */
export function obsColumnNamesForElement(
    spatialData: SpatialDataAssociationSource | undefined,
    elementType: AssociableSpatialElementType,
    elementKey: string,
): Set<string> | undefined {
    if (!spatialData) return undefined;
    const tables = spatialData.getAssociatedTables(spatialDataAssociationKind(elementType), elementKey);
    const table = tables.length === 1 ? tables[0][1] : undefined;
    const names = table?.getObsColumnNames?.();
    return names ? new Set(names) : undefined;
}

/**
 * Keep a column the viewer cannot read out of the viewer's inputs.
 *
 * Only for columns absent from the table's obs. Left in, the viewer would try to
 * load a column that is not there — a failed resolution and a user-visible notice
 * — and then be overruled by MDV's `featureState` colours anyway. Columns that ARE
 * in obs are handed over deliberately; see `obsColumnNamesForElement`.
 */
function withoutViewerFillColorColumn<T extends LayerConfig>(layer: T): T {
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
 *
 * Only the per-feature path needs this. A column the viewer owns is covered by the
 * viewer's own last-good retention, which holds the previous column's colours until
 * the new column's rows land; this covers MDV's asynchronous `loadColumnSet`, which
 * the viewer knows nothing about.
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
    const obsColumnNamesByAssociationKey = useMemo(() => {
        const entries: Array<[string, Set<string> | undefined]> = [];
        for (const elementKey of shapeKeys) {
            entries.push([`shapes:${elementKey}`, obsColumnNamesForElement(spatialData, "shapes", elementKey)]);
        }
        for (const elementKey of labelKeys) {
            entries.push([`labels:${elementKey}`, obsColumnNamesForElement(spatialData, "labels", elementKey)]);
        }
        return Object.fromEntries(entries);
    }, [spatialData, shapeKeys, labelKeys]);
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

    const layers = useMemo(() => {
        let changed = false;
        const nextLayers = { ...layerInputs.layers };

        for (const layerId of layerInputs.layerOrder) {
            const layer = layerInputs.layers[layerId];
            if (!isFillColorAssociableLayer(layer)) continue;

            const associationKey = `${layer.type}:${layer.elementKey}`;
            const table = tableByAssociationKey[associationKey];
            if (table?.status !== "resolved") continue;

            const fillColumnName = layerFillColorColumnName(layer);
            // Two routes, decided by whether the viewer can read the column at all.
            // In obs: hand it over with MDV's palette, and let the viewer do the
            // reading, the encoding and the load-window retention. Not in obs (a
            // gene score, `mdv_cell_id`, anything computed here): MDV is the only
            // party that has the values, so it colours each feature itself.
            const viewerFillColor =
                fillColumnName && obsColumnNamesByAssociationKey[associationKey]?.has(fillColumnName)
                    ? fillColorSchemeFromDataStore(table.dataStore, fillColumnName)
                    : undefined;
            const colorForRow =
                fillColumnName && !viewerFillColor
                    ? colorFunctionByDataSource[table.dataSourceName]?.[fillColumnName]
                    : undefined;
            // A selected column whose data has not loaded yet has no colour function.
            const colorReady = !fillColumnName || !!viewerFillColor || colorForRow !== undefined;
            const visibleRows =
                visibleRowsByDataSource[table.dataSourceName] ?? visibleRowsForDataStore(table.dataStore);

            // Filtering rides `featureState` on both routes: cross-filtering is
            // MDV's, not the table's, and no column in obs can express it.
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

            if (!viewerFillColor) {
                featureState = withPreservedFillColorsWhileLoading({
                    featureState,
                    fillColumnName,
                    colorReady,
                    previousFillColorByFeatureId: lastFillColorsRef.current[layerId],
                });
            }

            if (!fillColumnName || viewerFillColor) {
                // Nothing to preserve: either no column, or the viewer is holding
                // the last-good colours itself. Dropping the entry also stops a
                // previous MDV-drawn column bleeding into a later viewer-drawn one.
                const { [layerId]: _removed, ...remaining } = lastFillColorsRef.current;
                lastFillColorsRef.current = remaining;
            } else if (colorReady && featureState?.fillColorByFeatureId) {
                lastFillColorsRef.current[layerId] = featureState.fillColorByFeatureId;
            }

            if (!featureState && !viewerFillColor) continue;

            const projected = viewerFillColor
                ? { ...layer, fillColorByColumn: viewerFillColor }
                : withoutViewerFillColorColumn(layer);
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
        obsColumnNamesByAssociationKey,
        visibleRowsByDataSource,
        colorFunctionByDataSource,
    ]);

    return useMemo(
        () => ({
            layers,
            layerOrder: layerInputs.layerOrder,
        }),
        [layers, layerInputs.layerOrder],
    );
}
