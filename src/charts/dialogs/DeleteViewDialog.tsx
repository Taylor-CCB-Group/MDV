import {
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Typography,
} from "@mui/material";

const DeleteViewDialogComponent = (props: {
    open: boolean;
    onClose: () => void;
}) => {
    const cm = window.mdv.chartManager;

    const {
        viewManager,
        getState,
        _callListeners,
        changeView,
        removeAllCharts,
        showAddViewDialog,
    } = cm;
    let { viewData } = cm;

    const onDelete = () => {
        //todo - refactor this function to eliminate the usage of apply function
        //remove the view choice and change view to the next one
        const view = viewManager.current_view;

        // update the views
        const updatedViews = viewManager.all_views.filter((v) => v !== view);
        viewManager.setAllViews(updatedViews);

        const state = getState.apply(cm, []);

        //want to delete view and update any listeners
        state.view = null;

        _callListeners.apply(cm, ["state_saved", state]);

        if (updatedViews.length > 0) {
            // set current view to initial view
            const nextView = updatedViews[0];
            viewManager.setView(nextView);
            changeView.apply(cm, [nextView]);
        } else {
            // no other views exist
            removeAllCharts.apply(cm, []);
            viewData = {};
            showAddViewDialog.apply(cm, []);
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
            <DialogTitle>Delete View</DialogTitle>
            <DialogContent dividers>
                <Typography>Do you want to delete the current view?</Typography>
            </DialogContent>
            <DialogActions>
                <Button color="primary" onClick={onDelete}>
                    Yes
                </Button>
                <Button color="error" onClick={props.onClose}>
                    No
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default DeleteViewDialogComponent;
