import { useMemo } from "react";
import type { LoadedDataColumn, NumberDataType } from "@/charts/charts";
import { isDensityGridShown } from "../components/densityGridUtils";
import { type DualContourLegacyConfig, isLoadedNumericColumn } from "../contour_state";
import { useConfig, useFieldSpecs } from "../hooks";
import type { ChartArrayCell } from "./useChartArrayGrid";

export type DensityGridCells = {
    cells: ChartArrayCell[];
    densityFields: LoadedDataColumn<NumberDataType>[];
};

/**
 * One grid cell per density field that has loaded as a `double` column; a field still loading,
 * or of another type, has no cell. Loads nothing unless the chart is in grid mode.
 */
export function useDensityGridCells(): DensityGridCells {
    const config = useConfig<DualContourLegacyConfig>();
    const loadedFields = useFieldSpecs(isDensityGridShown(config) ? config.densityFields : undefined);
    const densityFields = useMemo(() => loadedFields.filter(isLoadedNumericColumn), [loadedFields]);
    const cells: ChartArrayCell[] = useMemo(
        () => densityFields.map((field) => ({ key: field.field, label: field.name || field.field })),
        [densityFields],
    );
    return { cells, densityFields };
}
