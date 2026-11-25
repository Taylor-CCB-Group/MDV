import GridViewOutlinedIcon from "@mui/icons-material/GridViewOutlined";
import { FormControlLabel, Switch, ToggleButton, ToggleButtonGroup } from "@mui/material";
import type { GridSplit } from "@/charts/charts";

/**
 * Under a column picker whose property the chart can be drawn as a grid of: whether it is, and if
 * so, a cell per column or a cell per category of the first. The picker takes several columns
 * only while it is on.
 */
export default function GridByControl({
    split,
    onChange,
}: {
    split: GridSplit | null;
    onChange: (split: GridSplit | null) => void;
}) {
    return (
        <div className="flex flex-wrap items-center gap-2">
            <FormControlLabel
                control={
                    <Switch
                        size="small"
                        checked={split !== null}
                        onChange={(_, checked) => onChange(checked ? "columns" : null)}
                    />
                }
                label="Grid by this"
                slotProps={{ typography: { variant: "caption" } }}
            />
            {split !== null && (
                <ToggleButtonGroup
                    size="small"
                    exclusive
                    value={split}
                    aria-label="Grid cells"
                    onChange={(_, next: GridSplit | null) => {
                        if (next) onChange(next);
                    }}
                >
                    <ToggleButton value="columns" sx={{ py: 0, textTransform: "none" }}>
                        A cell per column
                    </ToggleButton>
                    <ToggleButton value="categories" sx={{ py: 0, textTransform: "none" }}>
                        Per category
                    </ToggleButton>
                </ToggleButtonGroup>
            )}
        </div>
    );
}

/**
 * Marks the thing the chart is a grid of (in the layer list, the layer whose fill colour it is
 * split by), so it can be found without opening each layer's settings.
 */
export function GridByIcon({ split }: { split: GridSplit }) {
    const description =
        split === "columns"
            ? "The chart is a grid of this layer's colour: a cell per colour column"
            : "The chart is a grid of this layer's colour: a cell per category";
    return (
        <span className="flex flex-shrink-0" title={description} role="img" aria-label={description}>
            <GridViewOutlinedIcon fontSize="small" color="primary" sx={{ fontSize: 16 }} />
        </span>
    );
}
