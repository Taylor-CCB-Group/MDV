import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";
import type { ChartColumnImpact, ColumnRemovalImpact } from "@/charts/columnRemovalUtils";
import { describeColumnImpactReason } from "@/charts/columnRemovalUtils";
import {
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
        <Stack spacing={1.5}>
            <Typography fontWeight={600}>{title}</Typography>
            {charts.map((chart) => {
                const reasons = [...new Set(chart.reasons.map(describeColumnImpactReason))];
                return (
                    <Stack
                        key={chart.chartId}
                        direction={{ xs: "column", sm: "row" }}
                        justifyContent="space-between"
                        spacing={1}
                        sx={{
                            border: "1px solid",
                            borderColor: "divider",
                            borderRadius: 1,
                            p: 1.5,
                        }}
                    >
                        <Stack spacing={0.5}>
                            <Typography fontWeight={500}>
                                {chart.chartTitle || chart.chartTypeLabel}
                            </Typography>
                            <Typography variant="body2" color="text.secondary">
                                {chart.chartTypeLabel}
                            </Typography>
                        </Stack>
                        <Stack direction="row" spacing={0.75} useFlexGap flexWrap="wrap">
                            {reasons.map((reason) => (
                                <Chip
                                    key={`${chart.chartId}-${reason}`}
                                    size="medium"
                                    label={reason}
                                />
                            ))}
                        </Stack>
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

    return (
        <Dialog open={open} onClose={onClose} fullWidth maxWidth="md">
            <DialogTitle>
                Remove column "{columnName}"?
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers>
                <Stack spacing={2}>
                    <Typography variant="body2" color="text.secondary">
                        This column is used by other charts. Removing it will update or delete the charts below.
                    </Typography>
                    <ImpactSection title="Charts that will be updated" charts={updatedCharts} />
                    {updatedCharts.length > 0 && deletedCharts.length > 0 ? <Divider /> : null}
                    <ImpactSection title="Charts that will be deleted" charts={deletedCharts} />
                </Stack>
            </DialogContent>
            <DialogActions>
                <Button onClick={onClose}>Cancel</Button>
                <Button color={deletedCharts.length > 0 ? "error" : "primary"} onClick={onConfirm}>
                    Remove Column
                </Button>
            </DialogActions>
        </Dialog>
    );
}
