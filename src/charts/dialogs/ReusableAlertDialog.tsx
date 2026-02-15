import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";
import {
    Button,
    Container,
    Dialog,
    DialogActions,
} from "@mui/material";

export interface ReusableAlertDialogProps {
    open: boolean;
    handleClose: () => void;
    component: JSX.Element;
    isAlertErrorComponent?: boolean;
    isConfirmButton?: boolean;
    confirmText?: string;
    onConfirmClick?: () => void;
}

const ReusableAlertDialog = ({
    open,
    handleClose,
    component,
    isAlertErrorComponent,
    isConfirmButton,
    confirmText,
    onConfirmClick,
}: ReusableAlertDialogProps) => {
    return (
        <Dialog
            open={open}
            onClose={handleClose}
            disableRestoreFocus
            fullWidth
            maxWidth={isAlertErrorComponent ? "sm" : "md"}
            PaperProps={{
                sx: {
                    padding: 0,
                    overflow: "hidden",
                },
            }}
        >
            <DialogCloseIconButton onClose={handleClose} />
            <Container sx={{ mt: 1 }}>{component}</Container>
            {isConfirmButton && (
                <DialogActions>
                    <Button onClick={() => onConfirmClick?.()}>{confirmText ?? "Confirm"}</Button>
                </DialogActions>
            )}
        </Dialog>
    );
};

export default ReusableAlertDialog;
