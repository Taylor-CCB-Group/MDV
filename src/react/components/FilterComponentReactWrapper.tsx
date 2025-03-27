import { Autocomplete, TextField } from "@mui/material";
import { observer } from "mobx-react-lite";
import { useEffect, useState } from "react";

export type DropdownType = {
    options: string[];
};

const FilterDropdown = observer(() => {
    const cm = window.mdv.chartManager;
    const { viewManager } = cm;
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
