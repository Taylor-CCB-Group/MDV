import type {
    ChartColumnImpact,
    ColumnRemovalImpact,
    SavedViewColumnImpact,
} from "@/charts/columnRemovalUtils";
import {
    Accordion,
    AccordionDetails,
    AccordionSummary,
    Alert,
    Box,
    Button,
    Divider,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Stack,
    Typography,
} from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";

function getUsagePathLabel(path: string) {
    if (path === "param") return "param";
    if (path === "color_by") return "color_by";
    if (path === "tooltip") return "tooltip";
    if (path === "background_filter") return "background_filter";
    return path;
}

function getFixInstructions(chart: ChartColumnImpact) {
    const fixes = chart.usagePaths.map((path) => {
        if (path === "param") {
            return "Remove this column from the chart data, or remove the chart if it is no longer needed.";
        }
        if (path === "color_by") {
            return "Open the chart settings and choose a different column for color-by.";
        }
        if (path === "tooltip") {
            return "Open the chart settings and remove this column from the tooltip fields.";
        }
        if (path === "background_filter") {
            return "Open the chart settings and clear or change the background filter column.";
        }
        return `Open the chart settings and update "${path}" so it no longer uses this column.`;
    });

    return [...new Set(fixes)];
}

function getUsageSummaryLabel(chart: ChartColumnImpact) {
    if (chart.usagePaths.length === 1) {
        return getUsagePathLabel(chart.usagePaths[0]);
    }
    return `${chart.usagePaths.length} places`;
}

function SectionHeader({ count, title }: { count: number; title: string }) {
    return (
        <Stack alignItems="center" direction="row" spacing={1.25}>
            <Typography fontSize="1.05rem" fontWeight={700}>
                {title}
            </Typography>
            <Box
                sx={{
                    minWidth: 32,
                    height: 32,
                    borderRadius: "999px",
                    bgcolor: "action.selected",
                    color: "text.primary",
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                    px: 1,
                    fontSize: "0.95rem",
                    fontWeight: 700,
                }}
            >
                {count}
            </Box>
        </Stack>
    );
}

function ChartImpactCard({
    chart,
    onOpenView,
    viewName,
}: {
    chart: ChartColumnImpact;
    onOpenView?: (viewName: string) => void;
    viewName?: string;
}) {

    return (
        <Stack
            spacing={1.1}
            sx={{
                border: "1px solid",
                borderColor: "action.selected",
                borderRadius: 1.25,
                px: 2,
                py: 1.75,
                bgcolor: (theme) =>
                    theme.palette.mode === "dark" ? theme.palette.grey[800] : theme.palette.grey[50],
            }}
        >
            <Typography fontWeight={700}>
                {chart.chartTitle}
            </Typography>
            <Typography color="text.secondary" variant="body1">
                {chart.chartTypeLabel}
            </Typography>
            <Accordion
                disableGutters
                elevation={0}
                sx={{
                    border: "1px solid",
                    borderColor: "divider",
                    borderRadius: 1,
                    bgcolor: "transparent",
                    backgroundColor: "transparent",
                    boxShadow: "none",
                    "&::before": {
                        display: "none",
                    },
                }}
            >
                <AccordionSummary
                    expandIcon={<ExpandMoreIcon />}
                    sx={{
                        minHeight: "36px",
                        px: 1,
                        borderRadius: 1,
                        bgcolor: "action.hover",
                        "&.Mui-expanded": {
                            minHeight: "36px",
                        },
                        "& .MuiAccordionSummary-content": {
                            my: 0.5,
                        },
                        "& .MuiAccordionSummary-content.Mui-expanded": {
                            my: 0.5,
                        },
                        "& .MuiAccordionSummary-expandIconWrapper": {
                            color: "text.secondary",
                        },
                    }}
                >
                    <Typography fontWeight={700}>
                        How to fix
                    </Typography>
                </AccordionSummary>
                <AccordionDetails
                    sx={{
                        px: 1,
                        pt: 1,
                        pb: 0.25,
                        bgcolor: "transparent",
                    }}
                >
                    <Stack spacing={1.25}>
                        <Stack spacing={0.35}>
                            <Typography fontWeight={500} variant="body2">
                                Used as
                            </Typography>
                            <Typography color="text.secondary" variant="body2">
                                {chart.usagePaths.map((path) => getUsagePathLabel(path)).join(", ")}
                            </Typography>
                        </Stack>
                        <Stack spacing={0.35}>
                            <Typography fontWeight={500} variant="body2">
                                To fix
                            </Typography>
                            {getFixInstructions(chart).map((instruction) => (
                                <Typography key={instruction} color="text.secondary" variant="body2">
                                    {instruction}
                                </Typography>
                            ))}
                            {viewName && onOpenView ? (
                                <Stack alignItems="center" direction="row" flexWrap="wrap" spacing={0.5}>
                                    <Typography color="text.secondary" variant="body2">
                                        Go to the
                                    </Typography>
                                    <Button
                                        onClick={() => onOpenView(viewName)}
                                        size="small"
                                        sx={{ minWidth: 0, textTransform: "none", textDecoration: "underline" }}
                                        variant="text"
                                    >
                                        {viewName}
                                    </Button>
                                    <Typography color="text.secondary" variant="body2">
                                        view and update this chart there.
                                    </Typography>
                                </Stack>
                            ) : null}
                        </Stack>
                        <Stack spacing={0.35}>
                            <Typography fontWeight={500} variant="body2">
                                Then
                            </Typography>
                            {viewName && onOpenView ? (
                                <Typography color="text.secondary" variant="body2">
                                    After fixing the charts, save the view and try deleting the column again.
                                </Typography>
                            ) : (
                                <Typography color="text.secondary" variant="body2">
                                    After fixing the charts in the current view, try deleting the column again.
                                </Typography>
                            )}
                        </Stack>
                    </Stack>
                </AccordionDetails>
            </Accordion>
        </Stack>
    );
}

function ChartImpactList({
    charts,
    title,
}: {
    charts: ChartColumnImpact[];
    title: string;
}) {
    if (charts.length === 0) {
        return null;
    }

    return (
        <Stack spacing={1}>
            <SectionHeader count={charts.length} title={title} />
            {charts.map((chart) => (
                <ChartImpactCard
                    key={`${chart.chartId ?? chart.chartTitle}-${chart.chartType}-${chart.action}`}
                    chart={chart}
                />
            ))}
        </Stack>
    );
}

function SavedViewImpactSection({
    onOpenView,
    savedView,
}: {
    onOpenView: (viewName: string) => void;
    savedView: SavedViewColumnImpact;
}) {
    return (
        <Box
            sx={{
                border: "1px solid",
                borderColor: "divider",
                borderRadius: 1.25,
                px: 1.75,
                py: 1.5,
            }}
        >
            <Stack spacing={1.25}>
                <Stack
                    alignItems={{ sm: "center", xs: "flex-start" }}
                    direction={{ sm: "row", xs: "column" }}
                    justifyContent="space-between"
                    spacing={1}
                >
                    <Typography fontWeight={700}>{savedView.viewName}</Typography>
                    <Box
                        sx={{
                            minWidth: 32,
                            height: 32,
                            borderRadius: "999px",
                            bgcolor: "action.selected",
                            color: "text.primary",
                            display: "flex",
                            alignItems: "center",
                            justifyContent: "center",
                            px: 1,
                            fontSize: "0.95rem",
                            fontWeight: 700,
                        }}
                    >
                        {savedView.charts.length}
                    </Box>
                </Stack>
                <Divider />
                <Stack spacing={1}>
                    {savedView.charts.map((chart) => (
                        <ChartImpactCard
                            key={`${savedView.viewName}-${chart.chartId ?? chart.chartTitle}-${chart.chartType}-${chart.action}`}
                            chart={chart}
                            onOpenView={onOpenView}
                            viewName={savedView.viewName}
                        />
                    ))}
                </Stack>
            </Stack>
        </Box>
    );
}

function SavedViewImpactGroup({
    onOpenView,
    savedViews,
    title,
}: {
    onOpenView: (viewName: string) => void;
    savedViews: SavedViewColumnImpact[];
    title: string;
}) {
    if (savedViews.length === 0) {
        return null;
    }

    return (
        <Stack spacing={1.5}>
            <SectionHeader count={savedViews.length} title={title} />
            {savedViews.map((savedView) => (
                <SavedViewImpactSection
                    key={`${title}-${savedView.viewName}`}
                    onOpenView={onOpenView}
                    savedView={savedView}
                />
            ))}
        </Stack>
    );
}

export default function ColumnRemovalImpactDialog({
    impact,
    onClose,
    onConfirm,
    onOpenView,
    open,
}: {
    impact: ColumnRemovalImpact | null;
    onClose: () => void;
    onConfirm: () => void;
    onOpenView: (viewName: string) => void;
    open: boolean;
}) {
    const currentViewCharts = impact?.currentViewCharts ?? [];
    const savedViews = impact?.savedViews ?? [];
    const hasAnyImpacts = currentViewCharts.length > 0 || savedViews.length > 0;
    const isBlocked = hasAnyImpacts;

    return (
        <Dialog
            fullWidth
            maxWidth={isBlocked ? "lg" : "md"}
            onClose={onClose}
            open={open}
            PaperProps={{
                sx: {
                    width: "min(1120px, calc(100vw - 32px))",
                },
            }}
        >
            <DialogTitle>
                {impact ? `Delete Column "${impact.columnName}"` : "Delete Column"}
            </DialogTitle>
            <DialogContent dividers>
                <Stack spacing={2}>
                    <Alert severity={"error"}>
                        <Typography fontWeight={700} variant="body1">
                            {isBlocked
                                ? "This column is currently used by charts or saved views."
                                : "This action is irreversible."}
                        </Typography>
                        <Typography sx={{ mt: 0.75 }} variant="body1">
                            {isBlocked
                                ? "Deletion is blocked for now. For each chart below, check the 'Used as' section to see exactly whether this column is used in param, color_by, tooltip, background_filter, or another chart setting."
                                : "Deleting this column will permanently remove it from the datasource. This action cannot be undone, and the current view will be saved immediately after deletion."}
                        </Typography>
                    </Alert>

                    {isBlocked && (
                        <Box
                            sx={{
                                border: "1px solid",
                                borderColor: "divider",
                                borderRadius: 1.5,
                                px: 1.75,
                                py: 1.5,
                            }}
                        >
                            <Stack spacing={0.75}>
                                <Typography fontWeight={700} variant="body2">
                                    What you need to do
                                </Typography>
                                <Typography color="text.secondary" variant="body2">
                                    1. Find the affected view and chart below.
                                </Typography>
                                <Typography color="text.secondary" variant="body2">
                                    2. Read the 'Used as' section to see exactly where the column is used.
                                </Typography>
                                <Typography color="text.secondary" variant="body2">
                                    3. Open that chart, remove this column from the listed place, then try deleting the column again.
                                </Typography>
                            </Stack>
                        </Box>
                    )}

                    <ChartImpactList charts={currentViewCharts} title="Current View" />

                    <SavedViewImpactGroup
                        onOpenView={onOpenView}
                        savedViews={savedViews}
                        title="Other Views"
                    />
                </Stack>
            </DialogContent>
            <DialogActions>
                <Button onClick={onClose}>Cancel</Button>
                <Button color="error" disabled={isBlocked} onClick={onConfirm}>
                    {isBlocked ? "Deletion Blocked" : "Permanently Delete Column"}
                </Button>
            </DialogActions>
        </Dialog>
    );
}
