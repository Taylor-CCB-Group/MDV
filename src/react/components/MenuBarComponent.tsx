import {
    Add as AddIcon,
    CloudUpload as CloudUploadIcon,
    Home as HomeIcon,
    Remove as RemoveIcon,
    SaveAs as SaveAsIcon,
    Save as SaveIcon,
} from "@mui/icons-material";
import { AppBar, Box, Chip, Divider, Toolbar, Typography } from "@mui/material";
import { LockIcon, LockOpenIcon } from "lucide-react";
import { observer } from "mobx-react-lite";
import ToggleThemeWrapper from "@/charts/dialogs/ToggleTheme";
import ViewDialogWrapper from "@/charts/dialogs/ViewDialogWrapper";
import { useCurrentProject, useProject } from "@/modules/ProjectContext";
import { buildDashboardUrl, getApiRootFromDir } from "@/utils/mdvRouting";
import { useChartManager } from "../hooks";
import ChatButtons from "./ChatButtons";
import CustomTooltip from "./CustomTooltip";
import DebugButton from "./DebugButton";
import IconWithTooltip from "./IconWithTooltip";
import ViewSelectorWrapper from "./ViewSelectorComponent";
import ViewThumbnailComponent from "./ViewThumbnailComponent";

const MenuBarComponent = observer(() => {
    const { viewManager, config } = useChartManager();
    const { mainApiRoute, root } = useProject();
    const projectDisplayName = useCurrentProject()?.name;

    const isEditable = config.permission === "edit";
    const ariaLabel = isEditable ? "editable" : "view only";
    const homeHref = buildDashboardUrl(mainApiRoute === "/" ? getApiRootFromDir(root) : mainApiRoute);

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

    const handleAddDataSourceClick = async () => {
        const { default: FileUploadDialogReact } = await import("@/charts/dialogs/FileUploadDialogWrapper");
        new FileUploadDialogReact();
    };

    return (
        <>
            <AppBar position="static" color="default" sx={{ boxShadow: "none" }}>
                <Toolbar className="ciview-main-menu-bar" sx={{ padding: 2 }} disableGutters>
                    {/* the left group gives way first, so the right-hand actions stay on screen */}
                    <Box sx={{ display: "flex", alignItems: "center", justifyContent: "flex-start", minWidth: 0 }}>
                        <IconWithTooltip
                            tooltipText="Back to Catalog"
                            href={homeHref}
                            iconButtonProps={{ "aria-label": "Back to Catalog" }}
                        >
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

                        {projectDisplayName && (
                            <>
                                <Divider
                                    orientation="vertical"
                                    flexItem
                                    sx={{ mx: 1.5, my: 1, display: { xs: "none", md: "block" } }}
                                />
                                <CustomTooltip tooltipText={projectDisplayName}>
                                    <Typography
                                        variant="subtitle1"
                                        component="span"
                                        noWrap
                                        aria-label={`Project: ${projectDisplayName}`}
                                        sx={{
                                            display: { xs: "none", md: "block" },
                                            fontSize: "0.9rem",
                                            fontWeight: 500,
                                            minWidth: 0,
                                            maxWidth: { md: "14rem", lg: "24rem" },
                                            px: 1.25,
                                            py: 0.25,
                                            border: 1,
                                            borderColor: "divider",
                                            borderRadius: 1,
                                        }}
                                    >
                                        {projectDisplayName}
                                    </Typography>
                                </CustomTooltip>
                            </>
                        )}
                    </Box>
                    <Box sx={{ display: "flex", alignItems: "center", flexShrink: 0 }}>
                        <CustomTooltip
                            tooltipText={
                                isEditable
                                    ? "You can edit views and add data."
                                    : "You can only view, editing is disabled."
                            }
                        >
                            <Box
                                role="img"
                                aria-label={ariaLabel}
                                tabIndex={0}
                                sx={{ marginRight: "10px", marginBottom: "2px" }}
                            >
                                {isEditable ? (
                                    <LockOpenIcon height={18} aria-hidden="true" focusable="false" />
                                ) : (
                                    <LockIcon height={18} aria-hidden="true" focusable="false" />
                                )}
                            </Box>
                        </CustomTooltip>
                        <ChatButtons />
                        <ToggleThemeWrapper />
                        <DebugButton />
                    </Box>
                </Toolbar>
            </AppBar>
        </>
    );
});

const MenuBarWrapper = () => {
    return <MenuBarComponent />;
};

export default MenuBarWrapper;
