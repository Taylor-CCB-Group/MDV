import type {
    ChartColumnImpact,
    ColumnRemovalImpact,
    SavedViewColumnImpact,
} from "@/charts/types/columnRemoval";
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
            return "Open the chart settings and remove or change the column from parameters, or remove the chart if it is no longer needed.";
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

function trimChartTitle(title: string, maxLength = 42) {
    if (title.length <= maxLength) {
        return title;
    }
    return `${title.slice(0, maxLength - 1)}…`;
}

function getChartHeader(chart: ChartColumnImpact) {
    const hasCustomTitle =
        chart.chartTitle.trim().length > 0 &&
        chart.chartTitle !== chart.chartType &&
        chart.chartTitle !== chart.chartTypeLabel;

    if (!hasCustomTitle) {
        return {
            primaryLabel: chart.chartTypeLabel,
            secondaryLabel: null,
            fullPrimaryLabel: chart.chartTypeLabel,
        };
    }

    return {
        primaryLabel: trimChartTitle(chart.chartTitle),
        secondaryLabel: chart.chartTypeLabel,
        fullPrimaryLabel: chart.chartTitle,
    };
}

function getUsageBadgeLabel(chart: ChartColumnImpact) {
    return `Used as ${chart.usagePaths.map((path) => getUsagePathLabel(path)).join(", ")}`;
}

function SectionHeader({ count, title }: { count: number; title: string }) {
    return (
        <Stack alignItems="center" direction="row" spacing={1.25}>
            <Typography fontWeight={700}>
                {title}
            </Typography>
            <Box
                sx={{
                    minWidth: 20,
                    height: 20,
                    borderRadius: "999px",
                    bgcolor: "action.selected",
                    color: "text.primary",
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                    px: 1,
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
    const { fullPrimaryLabel, primaryLabel, secondaryLabel } = getChartHeader(chart);

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
            <Stack
                alignItems={{ sm: "center", xs: "flex-start" }}
                direction={{ sm: "row", xs: "column" }}
                justifyContent="space-between"
                spacing={1}
            >
                <Stack minWidth={0} spacing={0.35} sx={{ flex: 1 }}>
                    <Stack alignItems="center" direction="row" minWidth={0} spacing={1}>
                        <Typography
                            fontWeight={700}
                            noWrap
                            title={fullPrimaryLabel}
                            sx={{ minWidth: 0 }}
                        >
                            {primaryLabel}
                        </Typography>
                        {secondaryLabel ? (
                            <Typography color="text.secondary" noWrap variant="body1">
                                ({secondaryLabel})
                            </Typography>
                        ) : null}
                    </Stack>
                </Stack>
                <Box
                    sx={{
                        borderRadius: "999px",
                        bgcolor: "action.hover",
                        px: 1,
                        py: 0.4,
                        fontSize: "0.8rem",
                        fontWeight: 500,
                        whiteSpace: "nowrap",
                    }}
                >
                    {getUsageBadgeLabel(chart)}
                </Box>
            </Stack>
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
                        p: 1,
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
                    <Typography fontWeight={700} variant="body2">
                        How to fix
                    </Typography>
                </AccordionSummary>
                <AccordionDetails
                    sx={{
                        px: 1,
                        pt: 1,
                        pb: 1,
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
                            minWidth: 20,
                            height: 20,
                            borderRadius: "999px",
                            bgcolor: "action.selected",
                            color: "text.primary",
                            display: "flex",
                            alignItems: "center",
                            justifyContent: "center",
                            px: 1,
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
                                ? "Deletion is blocked."
                                : "Are you sure you want to delete this column?"}
                        </Typography>
                        <Typography sx={{ mt: 0.75 }} variant="body1">
                            {isBlocked
                                ? "This column is still used by other charts or views. Use 'How to fix' below, then try deleting again."
                                : "Note: This is a soft delete. The column won't be visible in MDV, but it will still appear in exported datasource files."}
                        </Typography>
                    </Alert>

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
                    {isBlocked ? "Deletion Blocked" : "Delete Column & Save"}
                </Button>
            </DialogActions>
        </Dialog>
    );
}
