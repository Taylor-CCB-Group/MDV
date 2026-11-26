import GridViewOutlinedIcon from "@mui/icons-material/GridViewOutlined";
import { FormControlLabel, Switch } from "@mui/material";
import type { GridSplit } from "@/charts/charts";
import type { GriddedEntry } from "./layerFillGridUtils";

const CELLS: Record<GridSplit, string> = { columns: "a cell per column", categories: "a cell per category" };

/**
 * Under a column picker whose property the chart can be drawn as a grid of: whether it is. The
 * split follows the column's type (see `GridBySpec`), and the label says which it is, or would be:
 * a cell per category of a categorical column, or a cell per numeric column, when the picker
 * takes several.
 */
export default function GridByControl({
    split,
    available,
    onChange,
}: {
    /** How the chart is split now, or null when it is not a grid of this property. */
    split: GridSplit | null;
    /** How it would be split, or null when the column cannot be a grid. */
    available: GridSplit | null;
    onChange: (on: boolean) => void;
}) {
    const cells = split ?? available;
    return (
        <FormControlLabel
            control={
                <Switch
                    size="small"
                    checked={split !== null}
                    disabled={split === null && available === null}
                    onChange={(_, checked) => onChange(checked)}
                />
            }
            label={cells ? `Grid by this: ${CELLS[cells]}` : "Grid by this"}
            title={
                split === null && available === null
                    ? "Choose a numeric or categorical column to show a grid of it"
                    : undefined
            }
            slotProps={{ typography: { variant: "caption" } }}
        />
    );
}

/**
 * Marks what the chart is a grid of: in the layer list, the layer (or MDV overlay) the grid is made
 * from, or, given that layer's name, the collapsed group that holds it.
 */
export function GridByIcon({ gridded, layerName }: { gridded: GriddedEntry; layerName?: string }) {
    const layer = layerName ? `'${layerName}' (in this group)` : "this layer";
    const of = gridded.byColour ? `the colour of ${layer}` : layer;
    const description = `The chart is a grid of ${of}: ${gridded.cells}`;
    return (
        <span className="flex shrink-0" title={description} role="img" aria-label={description}>
            <GridViewOutlinedIcon fontSize="small" color="primary" sx={{ fontSize: 16 }} />
        </span>
    );
}
