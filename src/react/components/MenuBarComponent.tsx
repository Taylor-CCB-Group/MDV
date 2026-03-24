import { AppBar, Box, Chip, Toolbar } from "@mui/material";
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
import CustomTooltip from "./CustomTooltip";
import { LockIcon, LockOpenIcon } from "lucide-react";

const MenuBarComponent = () => {
    const [error, setError] = useState<DebugErrorComponentProps['error'] | null>(null);
    const [open, setOpen] = useState(false);
    const cm = window.mdv.chartManager;
    const { viewManager, config, containerDiv } = cm;
    const { root, mainApiRoute } = useProject();
    const { buildInfo } = useBuildInfo();

    const isEditable = config.permission === "edit";
    const ariaLabel = isEditable ? "editable" : "view only";

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
            const ve = window.mdv?.validationErrors;
            const hasValidationErrors =
                ve && ((ve.datasources?.length ?? 0) > 0 || (ve.charts?.length ?? 0) > 0);
            new DebugChartReactWrapper({
                chartTypes,
                datasources,
                views,
                state,
                buildInfo,
                ...(hasValidationErrors && { validationErrors: ve }),
            });
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
                        {isEditable && (
                            <>
                                <IconWithTooltip tooltipText="Save View" onClick={handleSaveButtonClick}>
                                    <SaveIcon />
                                </IconWithTooltip>
                                <IconWithTooltip tooltipText="Save View As..." onClick={handleSaveAsButtonClick}>
                                    <SaveAsIcon />
                                </IconWithTooltip>
                            </>
                        )}
                        {isEditable && config.all_views && (
                            <>
                                <IconWithTooltip tooltipText="Create New View" onClick={handleCreateViewClick}>
                                    <AddIcon />
                                </IconWithTooltip>
                                <IconWithTooltip tooltipText="Delete Current View" onClick={handleDeleteViewClick}>
                                    <RemoveIcon />
                                </IconWithTooltip>
                            </>
                        )}
                        {isEditable && (
                            <IconWithTooltip tooltipText="Add Datasource" onClick={handleAddDataSourceClick}>
                                <CloudUploadIcon />
                            </IconWithTooltip>
                        )}
                    </Box>
                    <Box sx={{ display: "flex", alignItems: "center" }}>
                        <CustomTooltip 
                            tooltipText={
                                isEditable ? 
                                "You can edit views and add data." : 
                                "You can only view, editing is disabled."
                            }
                        >
                            <Box
                                role="img"
                                aria-label={ariaLabel}
                                tabIndex={0}
                                sx={{ marginRight: "10px", marginBottom: "2px" }}
                            >
                                {
                                    isEditable ? 
                                        <LockOpenIcon height={18} aria-hidden="true" focusable="false" /> :
                                        <LockIcon height={18} aria-hidden="true" focusable="false" /> 
                                }
                            </Box>
                        </CustomTooltip>
                        <ChatButtons />
                        <ToggleThemeWrapper />
                        <IconWithTooltip tooltipText="Debug / Report issue" onClick={handleDebugButtonClick}>
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
