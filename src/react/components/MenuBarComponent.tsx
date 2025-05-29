import { AppBar, Box, Button, ThemeProvider, Toolbar, useTheme } from "@mui/material";
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
import { ProjectProvider, useProject } from "@/modules/ProjectContext";
import DebugChartReactWrapper from "./DebugJsonDialogReactWrapper";
import ViewDialogWrapper from "@/charts/dialogs/ViewDialogWrapper";
import { useEffect, useRef, useState } from "react";
import ReusableDialog from "@/charts/dialogs/ReusableDialog";
import DebugErrorComponent, { type DebugErrorComponentProps } from "@/charts/dialogs/DebugErrorComponent";
import useBuildInfo from "@/catalog/hooks/useBuildInfo";
import ChatDialog from "@/charts/dialogs/ChatDialog";
import ChatLogDialog from "@/charts/dialogs/ChatLogDialog";
import ChatBubbleIcon from "@mui/icons-material/ChatBubble";
import ChatLogIcon from "@mui/icons-material/WebStories";
import { useChartManager } from "../hooks";
import { createPortal } from "react-dom";
export interface PopoutWindowProps {
    children: React.ReactNode;
    onClose: () => void;
    features?: string;
    title?: string;
  }

  export function PopoutWindow({
    children,
    features = "width=600,height=400",
    title = "",
    onClose,
  }: PopoutWindowProps) {
    const containerEl = useRef<HTMLElement>(document.createElement("div"));
    const externalWindow = useRef<Window | null>(null);

    useEffect(() => {
        externalWindow.current = window.open("", "_blank", features);
        if (!externalWindow.current) return;

        externalWindow.current.document.title = "Chat Window";
        // Function to synchronize the theme by syncing the class names
        const syncTheme = () => {
            if (!externalWindow.current) return;
            externalWindow.current.document.documentElement.className = document.documentElement.className;
        };

        // Initial sync
        syncTheme();

        // Set up a MutationObserver on the main window
        const themeObserver = new MutationObserver(syncTheme);
        themeObserver.observe(document.documentElement, {
            attributes: true,
            attributeFilter: ['class'],
        });

        // Function to add a stylesheet or style element to the popout window
        const addStyleElement = (
            styleElement: HTMLStyleElement | HTMLLinkElement,
        ) => {
            if (!externalWindow.current) {
                return;
            }
            if (styleElement instanceof HTMLLinkElement) {
                const newLink = externalWindow.current.document.createElement("link");
                newLink.rel = "stylesheet";
                newLink.href = styleElement.href;
                externalWindow.current.document.head.appendChild(newLink);
            } else if (styleElement instanceof HTMLStyleElement) {
                const newStyle = externalWindow.current.document.createElement("style");
                newStyle.innerHTML = styleElement.innerHTML;
                externalWindow.current.document.head.appendChild(newStyle);
            }
        };

        // Initial copy of all existing styles
        Array.from(
            document.querySelectorAll('link[rel="stylesheet"], style'),
        ).forEach((styleElement) => {
            addStyleElement(styleElement as HTMLStyleElement | HTMLLinkElement);
        });

        // Observe the main window's head for new stylesheets and style elements
        const observer = new MutationObserver((mutations) => {
            mutations.forEach((mutation) => {
                if (
                    mutation.type === "childList" &&
                    mutation.addedNodes.length > 0
                ) {
                    mutation.addedNodes.forEach((node) => {
                        if (
                            node instanceof HTMLLinkElement ||
                            node instanceof HTMLStyleElement
                        ) {
                            addStyleElement(node);
                        }
                    });
                }
            });
        });

        // Start observing the main window's head for changes
        observer.observe(document.head, {
            childList: true,
            subtree: true,
        });


        externalWindow.current.document.body.appendChild(containerEl.current);

        const handleUnload = () => {
            if (observer) observer.disconnect();
            if (themeObserver) themeObserver.disconnect();
            if (onClose) onClose();
        };
        externalWindow.current.addEventListener("beforeunload", handleUnload);

        // Cleanup
        return () => {
            if (externalWindow.current) {
                externalWindow.current.removeEventListener("beforeunload", handleUnload);
                externalWindow.current.close();
            }
        };
    }, [features, title, onClose]);

    return createPortal(children, containerEl.current);
}

const ChatButtons = () => {
    // very basic check to see if chat is enabled
    const chatEnabled = useChartManager().config.chat_enabled;
    const [open, setOpen] = useState(false);
    const [popout, setPopout] = useState(false);
    const theme = useTheme();
    const onClose = () => setOpen(false);
    const onPopout = () => {
        setOpen(false);
        setPopout(true);
    };

    const handleChatLogButtonClick = () => {
        new ChatLogDialog();
    };
    // if (!chatEnabled) return null;
    return (
        <>
            <IconWithTooltip tooltipText="Chat" onClick={() => setOpen(true)}>
                <ChatBubbleIcon />
            </IconWithTooltip>
            <IconWithTooltip tooltipText="Chat Log" onClick={handleChatLogButtonClick}>
                <ChatLogIcon />
            </IconWithTooltip>
            {!popout && <ChatDialog open={open} onClose={onClose} onPopout={onPopout} />}
            {
                popout && (
                    <PopoutWindow
                        onClose={() => {
                            onClose();
                            setPopout(false);
                        }}
                    >
                        <ProjectProvider>
                        <ThemeProvider theme={theme}>
                            <ChatDialog
                                open={true}
                                onClose={() => {
                                    onClose();
                                    setPopout(false);
                                }}
                                isPopout
                                fullscreen
                            />
                        </ThemeProvider>
                    </ProjectProvider>
                    </PopoutWindow>
                )
            }
        </>
    );
};

const MenuBarComponent = () => {
    const [error, setError] = useState<DebugErrorComponentProps['error'] | null>(null);
    const [open, setOpen] = useState(false);
    const cm = window.mdv.chartManager;
    const { viewManager, config, containerDiv } = cm;
    const { root } = useProject();
    const { buildInfo } = useBuildInfo();

    const handleHomeButtonClick = () => {
        // todo: uncomment when we fix state issues
        // viewManager.checkUnsavedState(() => {
            window.location.href = import.meta.env.DEV
                ? `${window.location.origin}/catalog_dev`
                : `${window.location.origin}/../`;
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
                <ReusableDialog
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
