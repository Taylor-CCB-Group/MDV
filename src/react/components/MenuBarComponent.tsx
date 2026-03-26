import { AppBar, Box, Chip, Toolbar } from "@mui/material";
import {
    Home as HomeIcon,
    Save as SaveIcon,
    SaveAs as SaveAsIcon,
    Add as AddIcon,
    Remove as RemoveIcon,
    CloudUpload as CloudUploadIcon,
} from "@mui/icons-material";
import ToggleThemeWrapper from "@/charts/dialogs/ToggleTheme";
import IconWithTooltip from "./IconWithTooltip";
import ViewSelectorWrapper from "./ViewSelectorComponent";
import ViewThumbnailComponent from "./ViewThumbnailComponent";
import FileUploadDialogReact from "@/charts/dialogs/FileUploadDialogWrapper";
import { useProject } from "@/modules/ProjectContext";
import ViewDialogWrapper from "@/charts/dialogs/ViewDialogWrapper";
import ChatButtons from "./ChatButtons";
import CustomTooltip from "./CustomTooltip";
import { LockIcon, LockOpenIcon } from "lucide-react";
import { observer } from "mobx-react-lite";
import { useChartManager } from "../hooks";
import DebugButton from "./DebugButton";

const MenuBarComponent = observer(() => {
    const { viewManager, config } = useChartManager();
    const { mainApiRoute } = useProject();

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
