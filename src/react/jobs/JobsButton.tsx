import { useProject } from "@/modules/ProjectContext";
import { Science as ScienceIcon } from "@mui/icons-material";
import { Box, Dialog, DialogContent, DialogTitle } from "@mui/material";
import { useState } from "react";
import IconWithTooltip from "../components/IconWithTooltip";
import { useChartManager, useDataSources } from "../hooks";
import JobsPanel from "./JobsPanel";
import { toJobsDataSource } from "./jobForm";
import { useJobColumnSync } from "./useJobColumnSync";

function JobsDialogContent() {
    const { root } = useProject();
    const sources = useDataSources().map(toJobsDataSource);
    return <JobsPanel root={root} sources={sources} />;
}

/** Menu bar button that opens the jobs selector (ADR-0012). */
export default function JobsButton() {
    const [open, setOpen] = useState(false);
    // the button stays mounted, so a job that finishes after the dialog closes still adds its columns
    useJobColumnSync(useProject().root, useChartManager());
    return (
        <>
            <IconWithTooltip tooltipText="Run Analysis" onClick={() => setOpen(true)}>
                <ScienceIcon />
            </IconWithTooltip>
            <Dialog open={open} onClose={() => setOpen(false)} fullWidth maxWidth="md">
                <DialogTitle>Run analysis</DialogTitle>
                <DialogContent>
                    {/* MUI drops DialogContent's top padding after a DialogTitle, which clips the first field's label */}
                    <Box sx={{ pt: 1 }}>{open && <JobsDialogContent />}</Box>
                </DialogContent>
            </Dialog>
        </>
    );
}
