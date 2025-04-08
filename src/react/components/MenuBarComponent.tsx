import { AppBar, Box, IconButton, Toolbar, Tooltip } from "@mui/material";
import {
    Home as HomeIcon,
    Save as SaveIcon,
    Add as AddIcon,
    Remove as RemoveIcon,
    CloudUpload as CloudUploadIcon,
    CameraAlt as CameraAltIcon,
    PestControl as PestControlIcon,
    BugReport as BugReportIcon,
} from "@mui/icons-material";
import ToggleThemeWrapper from "@/charts/dialogs/ToggleTheme";
import IconWithTooltip from "./IconWithTooltip";
import ViewSelectorWrapper from "./ViewSelectorComponent";
import ViewThumbnailComponent from "./ViewThumbnailComponent";
import FileUploadDialogReact from "@/charts/dialogs/FileUploadDialogWrapper";
import { toPng } from "html-to-image";
import { fetchJsonConfig } from "@/dataloaders/DataLoaderUtil";
import BaseChart from "@/charts/BaseChart";
import { getBuildInfo, useProject } from "@/modules/ProjectContext";
import DebugChartReactWrapper from "./DebugJsonDialogReactWrapper";
import ViewDialogWrapper from "@/charts/dialogs/ViewDialogWrapper";

const MenuBarComponent = () => {
    const cm = window.mdv.chartManager;
    const { viewManager, config, containerDiv } = cm;
    const { root } = useProject();

    const handleHomeButtonClick = () => {
        // const state = this.getState();
        // this._callListeners("state_saved", state);
        viewManager.checkUnsavedState(() => {
            window.location.href = import.meta.env.DEV
                ? `${window.location.origin}/catalog_dev`
                : `${window.location.origin}/../`;
        });
    };

    const handleSaveButtonClick = async () => {
        await viewManager.saveView();
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
        const buildInfo = getBuildInfo();
        new DebugChartReactWrapper({ chartTypes, datasources, views, state, buildInfo });
    };

    return (
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
                        <IconWithTooltip tooltipText="Save View" onClick={handleSaveButtonClick}>
                            <SaveIcon />
                        </IconWithTooltip>
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
                    <ToggleThemeWrapper />
                    <IconWithTooltip tooltipText="View Datasource Metadata" onClick={handleDebugButtonClick}>
                        <PestControlIcon sx={{height: "1.5rem", width: "1.5rem"}} />
                    </IconWithTooltip>
                </Box>
            </Toolbar>
        </AppBar>
    );
};

const MenuBarWrapper = () => {
    return <MenuBarComponent />;
};

export default MenuBarWrapper;
