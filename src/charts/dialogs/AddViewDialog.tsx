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
import { useEffect, useMemo, useState } from "react";

const AddViewDialogComponent = (props: {
    open: boolean;
    onClose: () => void;
}) => {
    const [viewName, setViewName] = useState("");
    const [isCloneView, setIsCloneView] = useState(false);
    const [checkedDs, setCheckedDs] = useState<{ [name: string]: boolean }>({});
    //must be a better way than using the window object
    const { viewManager, viewData } = window.mdv.chartManager;
    //we want a choice of all datasources (not just ones on the current view)
    //not sure how to get all the datasources available - this is just a temporary hack
    //if we don't use useMemo get constant re-rendering ?
    const dataSources = useMemo(()=>viewManager.cm.dataSources.map(x=>x.name),[viewManager.cm])
    const currentDataSources = useMemo(() => Object.keys(viewData.dataSources), [viewData.dataSources]);

    useEffect(() => {
        // Initialising the checkedDs with the data sources inside the view
        const tempDs: { [name: string]: boolean } = {};
        dataSources?.forEach?.((ds) => {
            tempDs[ds] =currentDataSources.indexOf(ds) !=-1;
        });
        setCheckedDs(tempDs);
    }, [dataSources]);

    const handleCheckboxChange = (ds: string) => {
        setCheckedDs((prevDs) => ({ ...prevDs, [ds]: !prevDs[ds] }));
    };
    
    const checkInvalidName = (name: string) => {
        return viewManager.all_views.includes(name);
    };

    const checkDsSelected = () => {
        if (dataSources.length < 2) return true; 
        const checkedDsArray = Object.values(checkedDs)?.filter?.((ds) => ds === true);
        return checkedDsArray.length > 0;
    };

    const onCreate = async () => {
        await viewManager.addView(viewName, checkedDs, isCloneView);
        props.onClose();
    };
    return (
        <Dialog open={props.open} onClose={props.onClose} fullWidth maxWidth="xs">
            <DialogTitle>Add New View</DialogTitle>
            <DialogContent dividers>
                <FormControlLabel
                    control={<Checkbox checked={isCloneView} onChange={(e) => setIsCloneView(e.target.checked)} />}
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
                    helperText={checkInvalidName(viewName) && "Name already exists"}
                />

                {dataSources.length > 1 && (
                    <FormGroup sx={{ mt: 2 }}>
                        {dataSources?.map?.((ds) => (
                            <FormControlLabel
                                key={ds}
                                control={
                                    <Checkbox
                                        // This !! is to ensure that the component is controlled and only boolean value is accepted
                                        checked={!!checkedDs[ds]}
                                        onChange={() => handleCheckboxChange(ds)}
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
                <Button
                    color="primary"
                    onClick={onCreate}
                    disabled={checkInvalidName(viewName) || !viewName || (!checkDsSelected() && !isCloneView)}
                >
                    Create View
                </Button>
                <Button color="error" onClick={props.onClose}>
                    Cancel
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default AddViewDialogComponent;
