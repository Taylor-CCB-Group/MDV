import {
    Button,
    Checkbox,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    FormControlLabel,
    FormGroup,
    TextField,
} from "@mui/material";
import { useEffect, useState } from "react";

const AddViewDialogComponent = (props: {
    open: boolean;
    onClose: () => void;
}) => {
    const [viewName, setViewName] = useState("");
    const [isCloneView, setIsCloneView] = useState(false);
    const [checkedDs, setCheckedDs] = useState<{ [name: string]: boolean }>({});
    const cm = window.mdv.chartManager;
    const { viewManager, viewData, dsIndex, contentDiv } = cm;

    useEffect(() => {
        // Initialising the checkedDs with the data sources inside the view
        const tempDs: { [name: string]: boolean } = {};
        Object.keys(viewData.dataSources)?.forEach?.((ds) => {
            tempDs[ds] = false;
        });
        setCheckedDs(tempDs);
    }, [viewData.dataSources]);

    const handleCheckboxChange = (ds: string) => {
        setCheckedDs((prevDs) => ({ ...prevDs, [ds]: !prevDs[ds] }));
    };

    const checkInvalidName = (name: string) => {
        return viewManager.all_views.includes(name);
    };

    const onCreate = () => {
        viewManager.setAllViews([...viewManager.all_views, viewName]);

        // Optionally make it the current view
        viewManager.setView(viewName);
        if (!isCloneView) {
            //remove all charts and links
            for (const ds in viewData.dataSources) {
                if (viewData.dataSources[ds].layout === "gridstack") {
                    const d = dsIndex[ds];
                    if (!d) continue;
                    cm.gridStack.destroy(d);
                }
            }
            cm.removeAllCharts();
            viewData.links = [];
            const state = cm.getState();
            state.view.initialCharts = {};
            state.view.dataSources = {};
            //only one datasource
            if (Object.keys(viewData.dataSources)?.length === 1) {
                const name = Object.keys(viewData.dataSources)?.[0];
                state.view.initialCharts[name] = [];
                state.view.dataSources[name] = {};
            } else {
                for (const ds in dsIndex) {
                    if (checkedDs[ds]) {
                        state.view.initialCharts[ds] = [];
                        state.view.dataSources[ds] = {};
                    }
                }
            }
            cm._callListeners("state_saved", state);
            contentDiv.innerHTML = "";
            cm._init(state.view);
        } else {
            const state = cm.getState();
            console.log("state add new: ", state);
            cm._callListeners("state_saved", state);
            cm._init(state.view);
        }
        props.onClose();
    };
    return (
        <Dialog
            open={props.open}
            onClose={props.onClose}
            fullWidth
            maxWidth="xs"
        >
            <DialogTitle>Add New View</DialogTitle>
            <DialogContent dividers>
                <FormControlLabel
                    control={
                        <Checkbox
                            checked={isCloneView}
                            onChange={(e) => setIsCloneView(e.target.checked)}
                        />
                    }
                    label="Clone current view"
                />

                <TextField
                    label="Name"
                    variant="outlined"
                    fullWidth
                    margin="normal"
                    value={viewName}
                    onChange={(e) => setViewName(e.target.value)}
                    error={checkInvalidName(viewName)}
                    helperText={
                        checkInvalidName(viewName) && "Name already exists"
                    }
                />

                {Object.keys(viewData.dataSources).length > 1 && (
                    <FormGroup sx={{ mt: 2 }}>
                        {Object.keys(viewData.dataSources)?.map?.((ds) => (
                            <FormControlLabel
                                key={ds}
                                control={
                                    <Checkbox
                                        // This !! is to ensure that the component is controlled and only boolean value is accepted
                                        checked={!!checkedDs[ds]}
                                        onChange={() =>
                                            handleCheckboxChange(ds)
                                        }
                                        disabled={isCloneView}
                                    />
                                }
                                label={`Include ${ds}`}
                            />
                        ))}
                    </FormGroup>
                )}
            </DialogContent>
            <DialogActions>
                <Button onClick={props.onClose}>Cancel</Button>
                <Button
                    color="primary"
                    onClick={onCreate}
                    disabled={checkInvalidName(viewName) || !viewName}
                >
                    Create New View
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default AddViewDialogComponent;
