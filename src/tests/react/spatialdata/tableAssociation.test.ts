import {
    type DataSourceAssociationCandidate,
    buildAssociatedFeatureStateFromRowMap,
    buildAssociatedShapesFeatureState,
    getShapesTableAssociation,
    resolveAssociatedElementTable,
    withPreservedFillColorsWhileLoading,
    withoutViewerFillColorColumn,
} from "@/react/spatialdata/table_association";
import type { ShapesRenderData } from "@spatialdata/core";
import type { LayerConfig } from "@spatialdata/vis";
import { describe, expect, test } from "vitest";

function shapesRenderData(featureIds: string[], rowIndexByFeatureIndex: number[]): ShapesRenderData {
    return {
        kind: "js-polygons",
        geometryKind: "polygon",
        elementKey: "cells",
        featureIds,
        polygons: [],
        rowIndexByFeatureIndex: Int32Array.from(rowIndexByFeatureIndex),
    };
}

describe("SpatialData table association", () => {
    test("reports resolved shape associations when features map to current rows", () => {
        const association = getShapesTableAssociation(
            shapesRenderData(["a", "b", "c"], [0, -1, 2]),
            3,
            "cells",
            "cell_datasource",
        );

        expect(association).toEqual({
            status: "resolved",
            tableName: "cells",
            dataSourceName: "cell_datasource",
            matchedFeatureCount: 2,
            featureCount: 3,
        });
    });

    test("leaves unmatched shapes unassociated", () => {
        expect(
            getShapesTableAssociation(shapesRenderData(["a", "b"], [-1, 12]), 2, "cells", "cell_datasource"),
        ).toEqual({ status: "none" });
    });

    test("resolves an associated shapes element through SpatialData table provenance", () => {
        const spatialData = {
            getAssociatedTables: (kind: string): Array<[string, unknown]> =>
                kind === "shapes" ? [["cells", null]] : [],
        };
        const cellsDataStore = {
            name: "cells_by_region",
            config: {
                spatialdata_tables: {
                    tables: [{ table_name: "cells" }],
                },
            },
        };
        const dataSources: DataSourceAssociationCandidate[] = [
            {
                name: "grid",
                dataStore: { name: "grid", config: {} },
            },
            {
                name: "cells_by_region",
                dataStore: cellsDataStore,
            },
        ];

        expect(
            resolveAssociatedElementTable({
                spatialData,
                elementType: "shapes",
                elementKey: "cell_boundaries",
                dataSources,
            }),
        ).toMatchObject({
            status: "resolved",
            tableName: "cells",
            dataSourceName: "cells_by_region",
        });
    });

    test("resolves labels through the same element association path", () => {
        const spatialData = {
            getAssociatedTables: (kind: string): Array<[string, unknown]> =>
                kind === "labels" ? [["segmentation_table", null]] : [],
        };
        const dataSources: DataSourceAssociationCandidate[] = [
            {
                name: "segmentation",
                dataStore: {
                    config: {
                        spatialdata_tables: {
                            tables: [{ table_name: "segmentation_table" }],
                        },
                    },
                },
            },
        ];

        expect(
            resolveAssociatedElementTable({
                spatialData,
                elementType: "labels",
                elementKey: "cell_labels",
                dataSources,
            }),
        ).toMatchObject({
            status: "resolved",
            tableName: "segmentation_table",
            dataSourceName: "segmentation",
        });
    });

    test("resolves labels from MDV table_name column metadata", () => {
        const spatialData = {
            getAssociatedTables: (kind: string): Array<[string, unknown]> =>
                kind === "labels" ? [["cell_binned", null]] : [],
        };
        const dataSources: DataSourceAssociationCandidate[] = [
            {
                name: "cells",
                dataStore: {
                    config: {
                        columns: [
                            {
                                field: "region",
                                values: ["cell_labels"],
                            },
                            {
                                field: "table_name",
                                values: ["cell_binned"],
                            },
                        ],
                    },
                },
            },
        ];

        expect(
            resolveAssociatedElementTable({
                spatialData,
                elementType: "labels",
                elementKey: "cell_labels",
                dataSources,
            }),
        ).toMatchObject({
            status: "resolved",
            tableName: "cell_binned",
            dataSourceName: "cells",
        });
    });

    test("builds feature-id keyed colors and hidden feature ids from row state", () => {
        const featureState = buildAssociatedShapesFeatureState({
            renderData: shapesRenderData(["a", "b", "c"], [0, 1, -1]),
            visibleRows: Uint32Array.from([0]),
            rowCount: 2,
            alpha: 123,
            colorForRow: (rowIndex) => (rowIndex === 0 ? [10, 20, 30] : [40, 50, 60]),
        });

        expect(featureState).toEqual({
            fillColorByFeatureId: {
                a: [10, 20, 30, 123],
                b: [40, 50, 60, 123],
            },
            hiddenFeatureIds: ["b"],
        });
    });

    test("preserves explicit feature state while adding table filter state", () => {
        const featureState = buildAssociatedShapesFeatureState({
            renderData: shapesRenderData(["a", "b"], [0, 1]),
            visibleRows: Uint32Array.from([1]),
            rowCount: 2,
            alpha: 255,
            baseFeatureState: {
                hiddenFeatureIds: ["manual"],
                fillColorByFeatureId: {
                    manual: [1, 2, 3, 4],
                },
            },
        });

        expect(featureState).toEqual({
            fillColorByFeatureId: {
                manual: [1, 2, 3, 4],
            },
            hiddenFeatureIds: ["manual", "a"],
        });
    });

    test("builds label feature colours from a feature-id row map", () => {
        const featureState = buildAssociatedFeatureStateFromRowMap({
            rowIndexByFeatureId: new Map([
                ["1", 0],
                ["2", 1],
            ]),
            visibleRows: Uint32Array.from([0, 1]),
            rowCount: 2,
            alpha: 255,
            colorForRow: (rowIndex) => (rowIndex === 0 ? [9, 8, 7] : [6, 5, 4]),
        });

        expect(featureState).toEqual({
            fillColorByFeatureId: {
                "1": [9, 8, 7, 255],
                "2": [6, 5, 4, 255],
            },
        });
    });

    test("strips fillColorByColumn from viewer layer configs", () => {
        const layer = {
            type: "labels",
            id: "labels-a",
            elementKey: "cell_labels",
            visible: true,
            opacity: 1,
            fillColorByColumn: {
                columnName: "Leiden",
                mode: "categorical",
            },
            featureState: {
                fillColorByFeatureId: {
                    "1": [1, 2, 3, 255],
                },
            },
        } as LayerConfig;

        expect(withoutViewerFillColorColumn(layer)).toEqual({
            type: "labels",
            id: "labels-a",
            elementKey: "cell_labels",
            visible: true,
            opacity: 1,
            featureState: {
                fillColorByFeatureId: {
                    "1": [1, 2, 3, 255],
                },
            },
        });
        expect(withoutViewerFillColorColumn(layer)).not.toHaveProperty("fillColorByColumn");
    });

    test("keeps previous fill colours while a newly selected column is still loading", () => {
        const previous = {
            a: [10, 20, 30, 180] as [number, number, number, number],
        };

        expect(
            withPreservedFillColorsWhileLoading({
                featureState: { hiddenFeatureIds: ["b"] },
                fillColumnName: "Leiden",
                colorReady: false,
                previousFillColorByFeatureId: previous,
            }),
        ).toEqual({
            hiddenFeatureIds: ["b"],
            fillColorByFeatureId: previous,
        });

        expect(
            withPreservedFillColorsWhileLoading({
                featureState: {
                    fillColorByFeatureId: {
                        a: [1, 2, 3, 180],
                    },
                },
                fillColumnName: "Leiden",
                colorReady: true,
                previousFillColorByFeatureId: previous,
            }),
        ).toEqual({
            fillColorByFeatureId: {
                a: [1, 2, 3, 180],
            },
        });
    });
});
