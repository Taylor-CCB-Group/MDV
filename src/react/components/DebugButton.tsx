import { Box } from "@mui/material";
import { PestControl as PestControlIcon } from "@mui/icons-material";
import { useState } from "react";
import { observer } from "mobx-react-lite";
import IconWithTooltip from "./IconWithTooltip";
import { useProject } from "@/modules/ProjectContext";
import { fetchJsonConfig } from "@/dataloaders/DataLoaderUtil";
import BaseChart from "@/charts/BaseChart";
import DebugChartReactWrapper from "./DebugJsonDialogReactWrapper";
import ReusableAlertDialog from "@/charts/dialogs/ReusableAlertDialog";
import DebugErrorComponent, { type DebugErrorComponentProps } from "@/charts/dialogs/DebugErrorComponent";
import useBuildInfo from "@/catalog/hooks/useBuildInfo";
import { useChartManager } from "../hooks";

const DebugButton = observer(() => {
    const [error, setError] = useState<DebugErrorComponentProps["error"] | null>(null);
    const [open, setOpen] = useState(false);

    const { validationFindings } = useChartManager();
    const { root } = useProject();
    const { buildInfo } = useBuildInfo();

    const handleDebugButtonClick = async () => {
        setError(null);
        setOpen(false);
        try {
            const datasources = await fetchJsonConfig(`${root}/datasources.json`, root);
            const views = await fetchJsonConfig(`${root}/views.json`, root);
            const state = await fetchJsonConfig(`${root}/state.json`, root);

            const chartTypes = Object.entries(BaseChart.types).map(([k, v]) => {
                const { class: omit, ...props } = v;
                return [k, props];
            });

            new DebugChartReactWrapper({
                chartTypes,
                datasources,
                views,
                state,
                buildInfo,
                ...(validationFindings?.hasAny && { validationFindings: validationFindings.snapshot }),
            });
        } catch (err) {
            setError(err instanceof Error ? {
                message: err.message,
                stack: err?.stack,
            } : {
                message: "Error occurred while fetching JSON",
                stack: `${err}`,
            });
            setOpen(true);
            console.error("Error occurred while fetching JSON: ", err);
        }
    };

    const showDot = Boolean(validationFindings?.hasAny);

    return (
        <>
            <IconWithTooltip tooltipText="Debug / Report issue" onClick={handleDebugButtonClick}>
                <Box sx={{ position: "relative", display: "inline-flex", alignItems: "center" }}>
                    <PestControlIcon sx={{ height: "1.5rem", width: "1.5rem" }} />
                    {showDot && (
                        <Box
                            sx={{
                                position: "absolute",
                                top: 0,
                                right: 0,
                                width: 8,
                                height: 8,
                                borderRadius: "50%",
                                backgroundColor: "#f2c94c",
                                boxShadow: "0 0 0 2px rgba(0,0,0,0.15)",
                            }}
                        />
                    )}
                </Box>
            </IconWithTooltip>

            {error && (
                <ReusableAlertDialog
                    open={open}
                    handleClose={() => setOpen(false)}
                    component={<DebugErrorComponent error={error} />}
                    isAlertErrorComponent
                />
            )}
        </>
    );
});

export default DebugButton;

