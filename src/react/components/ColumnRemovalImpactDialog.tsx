import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";
import type { ChartColumnImpact, ColumnRemovalImpact } from "@/types/columnRemovalTypes";
import { describeColumnImpactReason } from "@/charts/columnRemovalUtils";
import {
    Alert,
    Button,
    Chip,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Divider,
    Stack,
    Typography,
} from "@mui/material";

type ImpactSectionProps = {
    title: string;
    charts: ChartColumnImpact[];
};

function ImpactSection({ title, charts }: ImpactSectionProps) {
    if (charts.length === 0) {
        return null;
    }

    return (
        <Stack spacing={1}>
            <Stack direction="row" spacing={1} alignItems="center">
                <Typography fontWeight={600}>{title}</Typography>
                <Chip size="small" label={charts.length} />
            </Stack>
            {charts.map((chart) => {
                const reasons = [...new Set(chart.reasons.map(describeColumnImpactReason))];
                return (
                    <Stack
                        key={chart.chartId}
                        spacing={0.5}
                        sx={{
                            border: "1px solid",
                            borderColor: "divider",
                            borderRadius: 1,
                            px: 1.5,
                            py: 1.25,
                        }}
                    >
                        <Typography fontWeight={500}>
                            {chart.chartTitle || chart.chartTypeLabel}
                        </Typography>
                        <Typography variant="body2" color="text.secondary">
                            {chart.chartTypeLabel}
                            {reasons.length > 0 ? ` · uses this column in ${reasons.join(", ")}` : ""}
                        </Typography>
                    </Stack>
                );
            })}
        </Stack>
    );
}

export type ColumnRemovalImpactDialogProps = {
    open: boolean;
    columnName: string;
    impact: ColumnRemovalImpact | null;
    onClose: () => void;
    onConfirm: () => void;
};

export default function ColumnRemovalImpactDialog({
    open,
    columnName,
    impact,
    onClose,
    onConfirm,
}: ColumnRemovalImpactDialogProps) {
    const visibleCharts = impact?.charts.filter((chart) => !chart.isSourceChart) ?? [];
    const updatedCharts = visibleCharts.filter((chart) => chart.action === "update");
    const deletedCharts = visibleCharts.filter((chart) => chart.action === "delete");
    const hasDependentCharts = visibleCharts.length > 0;
    const updateCount = updatedCharts.length;
    const deleteCount = deletedCharts.length;

    return (
        <Dialog open={open} onClose={onClose} fullWidth maxWidth="md">
            <DialogTitle>
                Remove column "{columnName}"?
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers>
                <Stack spacing={2}>
                    <Alert severity={deleteCount > 0 ? "error" : "warning"}>
                        <Typography fontWeight={600} variant="body2">
                            This is a project-wide change and will be saved immediately.
                        </Typography>
                        <Typography variant="body2">
                            {hasDependentCharts
                                ? deleteCount > 0
                                    ? `${deleteCount} chart${deleteCount === 1 ? "" : "s"} will be deleted and ${updateCount} chart${updateCount === 1 ? "" : "s"} will be updated.`
                                    : `${updateCount} chart${updateCount === 1 ? "" : "s"} will be updated.`
                                : "Other saved views that still reference this column will be repaired when they are loaded."}
                        </Typography>
                    </Alert>
                    <ImpactSection title="Charts that will be updated" charts={updatedCharts} />
                    {updatedCharts.length > 0 && deletedCharts.length > 0 ? <Divider /> : null}
                    <ImpactSection title="Charts that will be deleted" charts={deletedCharts} />
                </Stack>
            </DialogContent>
            <DialogActions>
                <Button onClick={onClose}>Cancel</Button>
                <Button color="error" onClick={onConfirm}>
                    Remove Column And Save
                </Button>
            </DialogActions>
        </Dialog>
    );
}
