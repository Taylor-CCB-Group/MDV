import { Close } from "@mui/icons-material";
import { Button, Container, Dialog, DialogActions, DialogContent, DialogTitle, IconButton, Paper } from "@mui/material";

export interface ReusableDialogProps {
    open: boolean;
    handleClose: () => void;
    component: React.JSX.Element;
    isAlertErrorComponent?: boolean;
    isConfirmButton?: boolean;
    confirmText?: string;
    onConfirmClick?: () => void;
}

const ReusableDialog = ({
    open,
    handleClose,
    component,
    isAlertErrorComponent,
    isConfirmButton,
    confirmText,
    onConfirmClick,
}: ReusableDialogProps) => {
    return (
        <Dialog open={open} onClose={handleClose} disableRestoreFocus fullWidth maxWidth={isAlertErrorComponent ? "sm" : "md"}>
            <DialogContent dividers>
                <div className="flex items-center justify-center">
                    <Container>{component}</Container>
                </div>
            </DialogContent>
            <DialogActions>
                {isConfirmButton && <Button onClick={() => onConfirmClick?.()}>{confirmText ?? "Confirm"}</Button>}
                <Button onClick={handleClose} color="error">
                    Close
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default ReusableDialog;
