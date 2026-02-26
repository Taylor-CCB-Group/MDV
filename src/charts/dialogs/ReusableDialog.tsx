import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";
import {
    Button,
    Container,
    Dialog,
    DialogActions,
    DialogContent,
} from "@mui/material";

export interface ReusableDialogProps {
    open: boolean;
    handleClose: () => void;
    component: JSX.Element;
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
        <Dialog open={open} onClose={handleClose} disableRestoreFocus fullWidth maxWidth={"md"}>
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
