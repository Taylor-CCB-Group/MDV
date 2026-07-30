import type { ShapesRenderData } from "@spatialdata/core";
import { describe, expect, test } from "vitest";

import {
    buildAssociatedShapesFeatureState,
    type DataSourceAssociationCandidate,
    getShapesTableAssociation,
    resolveAssociatedShapesTable,
} from "@/react/spatialdata/table_association";

function shapesRenderData(
    featureIds: string[],
    rowIndexByFeatureIndex: number[],
): ShapesRenderData {
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
            getShapesTableAssociation(
                shapesRenderData(["a", "b"], [-1, 12]),
                2,
                "cells",
                "cell_datasource",
            ),
        ).toEqual({ status: "none" });
    });

    test("resolves a shapes element through SpatialData table provenance", () => {
        const spatialData = {
            getAssociatedTables: (): Array<[string, unknown]> => [["cells", null]],
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
            resolveAssociatedShapesTable({
                spatialData,
                elementKey: "cell_boundaries",
                dataSources,
            }),
        ).toMatchObject({
            status: "resolved",
            tableName: "cells",
            dataSourceName: "cells_by_region",
        });
    });

    test("builds feature-id keyed colors and hidden feature ids from row state", () => {
        const featureState = buildAssociatedShapesFeatureState({
            renderData: shapesRenderData(["a", "b", "c"], [0, 1, -1]),
            visibleRows: Uint32Array.from([0]),
            rowCount: 2,
            alpha: 123,
            colorForRow: (rowIndex) =>
                rowIndex === 0 ? [10, 20, 30] : [40, 50, 60],
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
});
