import { useMemo } from "react";
import type { LoadedDataColumn, NumberDataType } from "@/charts/charts";
import { DensityOverlayLayer } from "@/webgl/SpatialLayer";
import { type DualContourLegacyConfig, useFieldContour } from "../contour_state";
import { getFieldColor } from "../fieldColorManager";
import { useConfig } from "../hooks";
import type { UseChartArrayGridResult } from "../hooks/useChartArrayGrid";
import { useRange } from "../spatial_context";
import type { DeckCellLabel } from "./DeckCellGrid";

/**
 * For a density grid: the cells in view, a density layer for each of them (its id prefixed with
 * the cell's view id, so it draws only there) and a label per cell.
 */
export function useDensityGridCellLayers(
    grid: UseChartArrayGridResult,
    densityFields: readonly LoadedDataColumn<NumberDataType>[],
) {
    const config = useConfig<DualContourLegacyConfig>();
    const { modelMatrix } = useRange();
    const visibleCells = useMemo(
        () => grid.visibleCellIndices.filter((index) => densityFields[index] && grid.viewIds[index]),
        [grid.visibleCellIndices, grid.viewIds, densityFields],
    );
    const visibleFields = useMemo(
        () => visibleCells.map((index) => densityFields[index]),
        [visibleCells, densityFields],
    );
    // aligned with `visibleFields`
    const contours = useFieldContour({
        id: "densityGrid",
        fields: visibleFields,
        fill: config.contour_fill,
        bandwidth: config.contour_bandwidth ?? 0.1,
        intensity: config.contour_intensity ?? 0.1,
        opacity: config.contour_opacity ?? 0.2,
        fillThreshold: config.contour_fillThreshold ?? 2,
    });
    const densityLayers = useMemo(
        () =>
            visibleCells.map(
                (index, i) =>
                    new DensityOverlayLayer({
                        id: `${grid.viewIds[index]}-density`,
                        contourLayers: [contours[i]],
                        modelMatrix,
                    }),
            ),
        [visibleCells, grid.viewIds, contours, modelMatrix],
    );
    const labels: DeckCellLabel[] = useMemo(
        () => densityFields.map((field) => ({ text: field.name || field.field, color: getFieldColor(field.field) })),
        [densityFields],
    );
    return { visibleCells, densityLayers, labels };
}
