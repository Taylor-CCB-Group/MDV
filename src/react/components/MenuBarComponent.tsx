import { AppBar, Box, Toolbar } from "@mui/material";
import {
    Home as HomeIcon,
    Save as SaveIcon,
    SaveAs as SaveAsIcon,
    Add as AddIcon,
    Remove as RemoveIcon,
    CloudUpload as CloudUploadIcon,
    PestControl as PestControlIcon,
} from "@mui/icons-material";
import ToggleThemeWrapper from "@/charts/dialogs/ToggleTheme";
import IconWithTooltip from "./IconWithTooltip";
import ViewSelectorWrapper from "./ViewSelectorComponent";
import ViewThumbnailComponent from "./ViewThumbnailComponent";
import FileUploadDialogReact from "@/charts/dialogs/FileUploadDialogWrapper";
import { fetchJsonConfig } from "@/dataloaders/DataLoaderUtil";
import BaseChart from "@/charts/BaseChart";
import { useProject } from "@/modules/ProjectContext";
import DebugChartReactWrapper from "./DebugJsonDialogReactWrapper";
import ViewDialogWrapper from "@/charts/dialogs/ViewDialogWrapper";
import { useState } from "react";
import DebugErrorComponent, { type DebugErrorComponentProps } from "@/charts/dialogs/DebugErrorComponent";
import useBuildInfo from "@/catalog/hooks/useBuildInfo";
import ChatButtons from "./ChatButtons";
import ReusableAlertDialog from "@/charts/dialogs/ReusableAlertDialog";

const MenuBarComponent = () => {
    const [error, setError] = useState<DebugErrorComponentProps['error'] | null>(null);
    const [open, setOpen] = useState(false);
    const cm = window.mdv.chartManager;
    const { viewManager, config, containerDiv } = cm;
    const { root, mainApiRoute } = useProject();
    const { buildInfo } = useBuildInfo();

    const handleHomeButtonClick = () => {
        // todo: uncomment when we fix state issues
        // viewManager.checkUnsavedState(() => {
            window.location.href = import.meta.env.DEV
                ? `${window.location.origin}/catalog_dev` //todo: sort out routing for dev
                : mainApiRoute;
        // });
    };

    const handleSaveButtonClick = async () => {
        await viewManager.saveView();
    };

    const handleSaveAsButtonClick = () => {
        new ViewDialogWrapper("save_as");
    };

    const handleCreateViewClick = async () => {
        viewManager.checkAndAddView();
    };

    const handleDeleteViewClick = () => {
        new ViewDialogWrapper("delete");
    };

    const handleAddDataSourceClick = () => {
        new FileUploadDialogReact();
    };

    const handleDebugButtonClick = async () => {
        setError(null);
        setOpen(false);
        try {
            const datasources = await fetchJsonConfig(
                `${root}/datasources.json`,
                root,
            );
            const views = await fetchJsonConfig(`${root}/views.json`, root);
            const state = await fetchJsonConfig(`${root}/state.json`, root);
            const chartTypes = Object.entries(BaseChart.types).map(([k, v]) => {
                const { class: omit, ...props } = v;
                return [k, props];
            });
            new DebugChartReactWrapper({ chartTypes, datasources, views, state, buildInfo });
        } catch (error) {
            setError(error instanceof Error ? {
                message: error.message,
                stack: error?.stack,
            } : {
                message: "Error occurred while fetching JSON",
                stack: `${error}`
            });
            setOpen(true);
            console.error("Error occurred while fetching JSON: ", error);
        }

    };

    return (
        <>
            <AppBar position="static" color="default" sx={{ boxShadow: "none" }}>
                <Toolbar className="ciview-main-menu-bar" sx={{ padding: 2 }} disableGutters>
                    <Box sx={{ display: "flex", alignItems: "center", justifyContent: "flex-start" }}>
                        <IconWithTooltip tooltipText="Back to Catalog" onClick={handleHomeButtonClick}>
                            <HomeIcon />
                        </IconWithTooltip>
                        {config.all_views && (
                            <Box>
                                <ViewSelectorWrapper />
                                <ViewThumbnailComponent />
                            </Box>
                        )}
                        {config.permission === "edit" && (
                            <>
                                <IconWithTooltip tooltipText="Save View" onClick={handleSaveButtonClick}>
                                    <SaveIcon />
                                </IconWithTooltip>
                                <IconWithTooltip tooltipText="Save View As..." onClick={handleSaveAsButtonClick}>
                                    <SaveAsIcon />
                                </IconWithTooltip>
                            </>
                        )}
                        {config.permission === "edit" && config.all_views && (
                            <>
                                <IconWithTooltip tooltipText="Create New View" onClick={handleCreateViewClick}>
                                    <AddIcon />
                                </IconWithTooltip>
                                <IconWithTooltip tooltipText="Delete Current View" onClick={handleDeleteViewClick}>
                                    <RemoveIcon />
                                </IconWithTooltip>
                            </>
                        )}
                        {config.permission === "edit" && (
                            <IconWithTooltip tooltipText="Add Datasource" onClick={handleAddDataSourceClick}>
                                <CloudUploadIcon />
                            </IconWithTooltip>
                        )}
                    </Box>
                    <Box sx={{ display: "flex", alignItems: "center" }}>
                        <ChatButtons />
                        <ToggleThemeWrapper />
                        <IconWithTooltip tooltipText="View Datasource Metadata" onClick={handleDebugButtonClick}>
                            <PestControlIcon sx={{height: "1.5rem", width: "1.5rem"}} />
                        </IconWithTooltip>
                    </Box>
                </Toolbar>
            </AppBar>
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
};

const MenuBarWrapper = () => {
    return <MenuBarComponent />;
};

export default MenuBarWrapper;
