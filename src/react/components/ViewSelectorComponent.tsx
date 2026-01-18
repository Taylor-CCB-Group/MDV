import { fetchJsonConfig } from "@/dataloaders/DataLoaderUtil";
import { useProject } from "@/modules/ProjectContext";
import { Autocomplete, IconButton, TextField } from "@mui/material";
import { PriorityHigh as PriorityHighIcon } from "@mui/icons-material";
import { observer } from "mobx-react-lite";
import { useEffect, useState } from "react";
import { useChartManager, useViewManager } from "../hooks";
import DebugErrorComponent, { type DebugErrorComponentProps } from "@/charts/dialogs/DebugErrorComponent";
import ReusableAlertDialog from "@/charts/dialogs/ReusableAlertDialog";

export type DropdownType = {
    options: string[];
};

function useUpdateViewList() {
    const viewManager = useViewManager();
    const { root } = useProject();

    useEffect(() => {
        const interval = setInterval(async () => {
            try {
                const config = await fetchJsonConfig(`${root}/state.json`, root);
                // we should consider zod validation here... but let's just try not to mess up
                if (!config.all_views) {
                    console.warn(`expected state.json, got ${config}`);
                    return;
                }
                viewManager.setAllViews(config.all_views);
            } catch (error) {
                //! not setting error here because this will be called every 0.5 sec which is causing the error dialog to appear multiple times
                console.error("Error fetching JSON: ", error);
            }
        }, 500);

        return () => {
            clearInterval(interval);
        };
    });
}

function useKeyboardShortcuts() {
    const viewManager = useViewManager();

    useEffect(() => {
        const handleKeyDown = (event: KeyboardEvent) => {
            if (event.key === "s" && (event.ctrlKey || event.metaKey)) {
                event.preventDefault();
                // Save the current view
                try {
                    viewManager.saveView();
                } catch (error) {
                    console.error("Error saving view:", error);
                    // we should consider an error dialog here
                }
            }
        };

        window.addEventListener("keydown", handleKeyDown);

        return () => {
            window.removeEventListener("keydown", handleKeyDown);
        };
    }, [viewManager]);
}

const ViewSelectorDropdown = observer(() => {
    const cm = useChartManager();
    const viewManager = useViewManager();

    useKeyboardShortcuts();
    // todo: uncomment when we fix state issues
    // useUpdateViewList();

    const options = viewManager.all_views;
    const [dirty, setDirty] = useState(false);
    const [error, setError] = useState<DebugErrorComponentProps['error']>();
    const [openError, setOpenError] = useState(false);

    // todo: uncomment when we fix state issues
    // useEffect(() => {
    //     const interval = setInterval(() => {
    //         setDirty((v) => {
    //             const dirty = viewManager.hasUnsavedChanges();
    //             if (!v && dirty) {
    //                 viewManager.hasUnsavedChanges(true);
    //             }
    //             return dirty;
    //         });
    //     }, 1000);

    //     return () => clearInterval(interval);
    // }, [viewManager]);

    // useEffect(() => {
    //     cm.addListener("view_selector", (type: string, data: any) => {
    //         if (type === "view_loaded") {
    //             setDirty(false);
    //         }
    //     });

    //     return () => {
    //         cm.removeListener("view_selector");
    //     };
    // }, [cm]);

    return (
        <>
            <Autocomplete
                options={options}
                value={viewManager.current_view || null}
                onChange={(_event, newValue) => {
                    if (newValue) {
                        cm.changeView(newValue);
                    }
                }}
                // todo: revert back when we fix state issues
                // renderInput={(params) => <TextField {...params} label={`Select View${dirty ? "*" : ""}`} />}
                renderInput={(params) => <TextField {...params} label={"Select View"} />}
                sx={{ display: "inline-flex", width: "20vw", margin: "0.2em" }}
            />
            {error && (
                <IconButton color="error" onClick={() => setOpenError(true)}>
                    <PriorityHighIcon />
                </IconButton>
            )}
            {openError && (
                <ReusableAlertDialog
                    open={openError}
                    handleClose={() => setOpenError(false)}
                    component={
                        <DebugErrorComponent
                            // todo: update the error to right format before passing
                            error={{ message: error?.message as string, stack: error?.stack }}
                            // extraMetadata={extraMetaData}
                        />
                    }
                />
            )}
        </>
    );
});

const ViewSelectorWrapper = () => {
    return <ViewSelectorDropdown />;
};

export default ViewSelectorWrapper;
