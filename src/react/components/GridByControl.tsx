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
