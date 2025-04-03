import { fetchJsonConfig } from "@/dataloaders/DataLoaderUtil";
import { useProject } from "@/modules/ProjectContext";
import { Autocomplete, TextField } from "@mui/material";
import { observer } from "mobx-react-lite";
import { useEffect, useState } from "react";
import { useChartManager, useViewManager } from "../hooks";

export type DropdownType = {
    options: string[];
};

function useUpdateViewList() {
    const viewManager = useViewManager();
    const { root } = useProject();
    useEffect(() => {
        const interval = setInterval(async () => {
            const config = await fetchJsonConfig(`${root}/state.json`, root);
            // we should consider zod validation here... but let's just try not to mess up
            if (!config.all_views) {
                console.warn(`expected state.json, got ${config}`);
                return;
            }
            viewManager.setAllViews(config.all_views);
        }, 500);
        return () => {
            clearInterval(interval);
        }
    });
}

function useKeyboardShortcuts() {
    const viewManager = useViewManager();
    useEffect(() => {
        const handleKeyDown = (event: KeyboardEvent) => {
            if (event.key === "s" && (event.ctrlKey || event.metaKey)) {
                event.preventDefault();
                // Save the current view
                viewManager.saveView();
            }
        };
        window.addEventListener("keydown", handleKeyDown);
        return () => {
            window.removeEventListener("keydown", handleKeyDown);
        };
    }, [viewManager]);
}


const FilterDropdown = observer(() => {
    const cm = useChartManager();
    const viewManager = useViewManager();
    useKeyboardShortcuts();
    useUpdateViewList();
    const options = viewManager.all_views;
    const [dirty, setDirty] = useState(false);
    useEffect(() => {
        const interval = setInterval(() => {
            setDirty(v => {
                const dirty = viewManager.hasUnsavedChanges();
                if (!v && dirty) {
                    viewManager.hasUnsavedChanges(true);
                }
                return dirty;
            });
        }, 1000);
        return () => clearInterval(interval);
    }, [viewManager]);
    useEffect(() => {
        cm.addListener("view_selector", (type: string, data: any) => {
            if (type === "view_loaded") {
                setDirty(false);
            }
        });
        return () => {
            cm.removeListener("view_selector");
        }
    }, [cm]);
    return (
        <Autocomplete
            options={options}
            value={viewManager.current_view || null}
            onChange={(_event, newValue) => {
                if (newValue) {
                    // Updating state and changing view
                    // const state = cm.getState();
                    // cm._callListeners("state_saved", state);
                    cm.changeView(newValue);
                }
            }}
            renderInput={(params) => (
                <TextField {...params} label={`Select View${dirty ? '*' : ''}`} />
            )}
            sx={{ display: "inline-flex", width: "20vw", margin: "0.2em" }}
        />
    );
});

const FilterComponentReactWrapper = () => {
    return <FilterDropdown />;
};

export default FilterComponentReactWrapper;
