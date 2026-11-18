import { useMemo } from "react";
import { type DualContourLegacyConfig, isLoadedDoubleColumn } from "../contour_state";
import { useConfig, useFieldSpecs } from "../hooks";
import type { ChartArrayCell } from "./useChartArrayGrid";

/** One grid cell per density field that has loaded; a field still loading has no cell yet. */
export function useDensityGridCells() {
    const config = useConfig<DualContourLegacyConfig>();
    const loadedFields = useFieldSpecs(config.densityFields);
    const densityFields = useMemo(() => loadedFields.filter(isLoadedDoubleColumn), [loadedFields]);
    const cells: ChartArrayCell[] = useMemo(
        () => densityFields.map((field) => ({ key: field.field, label: field.name || field.field })),
        [densityFields],
    );
    return {
        cells,
        densityFields,
        configuredFieldCount: config.densityFields?.length ?? 0,
    };
}
