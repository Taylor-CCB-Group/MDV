import { Close } from "@mui/icons-material";
import { Button, Container, Dialog, DialogActions, DialogContent, DialogTitle, IconButton, Paper } from "@mui/material";

export interface ReusableDialogProps {
    open: boolean;
    handleClose: () => void;
    component: JSX.Element;
    isErrorMessage?: boolean;
}

const ReusableDialog = ({
    open,
    handleClose,
    component,
    isErrorMessage,
}: ReusableDialogProps) => {
    return (
        <Dialog
            open={open}
            onClose={handleClose}
            disableRestoreFocus
            fullWidth
            maxWidth={isErrorMessage ? 'sm' : 'md'}
        >
            <DialogContent dividers>
            <div
                className="flex items-center justify-center"
            >
                    <Container>{component}</Container>
            </div>
            </DialogContent>
            <DialogActions>
                <Button onClick={handleClose}>Close</Button>
            </DialogActions>
        </Dialog>
    );
};

export default ReusableDialog;
