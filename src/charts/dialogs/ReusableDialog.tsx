import { Close } from "@mui/icons-material";
import { Container, Dialog, IconButton, Paper } from "@mui/material";

export interface ReusableDialogProps {
    open: boolean;
    handleClose: () => void;
    component: JSX.Element;
}

const ReusableDialog = ({
    open,
    handleClose,
    component,
}: ReusableDialogProps) => {
    return (
        <Dialog
            open={open}
            onClose={handleClose}
            fullScreen
            disableEscapeKeyDown={true}
            PaperProps={{
                style: {
                    backgroundColor: "var(--fade_background_color)",
                    backdropFilter: "blur(1px)",
                },
            }}
        >
            <IconButton
                onClick={handleClose}
                className="absolute top-4 right-4"
                sx={{
                    position: "absolute",
                    top: 16,
                    right: 16,
                    backgroundColor: "var(--background_color)",
                    "&:hover": {
                        backgroundColor: "var(--menu_bar_color)",
                    },
                }}
            >
                <Close />
            </IconButton>

            <div className="h-screen flex items-center justify-center">
                <Paper elevation={24} sx={{ p: 2 }}>
                    <Container>{component}</Container>
                </Paper>
            </div>
        </Dialog>
    );
};

export default ReusableDialog;
