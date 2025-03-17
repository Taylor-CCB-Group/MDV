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

    const { viewManager } = cm;
    let { viewData } = cm;

    const onDelete = () => {
        //remove the view choice and change view to the next one
        const view = viewManager.current_view;

        // update the views
        const updatedViews = viewManager.all_views.filter((v) => v !== view);
        viewManager.setAllViews(updatedViews);

        const state = cm.getState();

        //want to delete view and update any listeners
        state.view = null;

        cm._callListeners("state_saved", state);

        if (updatedViews.length > 0) {
            // set current view to initial view
            const nextView = updatedViews[0];
            viewManager.setView(nextView);
            cm.changeView(nextView);
        } else {
            // no other views exist
            cm.removeAllCharts();
            viewData = {};
            cm.showAddViewDialog();
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
