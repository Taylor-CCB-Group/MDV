import { Button, Dialog, DialogActions, DialogContent, DialogTitle, Link, Typography } from "@mui/material";
import { DialogCloseIconButton } from "./ProjectRenameModal";
import { HELP_DIALOG_MESSAGE, MDV_DOC_WEBSITE, MDV_EMAIL } from "@/utilities/constants";

export type HelpDialogType = {
    open: boolean;
    onClose: () => void;
};

const HelpDialog = ({ open, onClose }: HelpDialogType) => {
    return (
        <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
            <DialogTitle>
                Help/Feedback
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers>
                <Typography variant="body1" sx={{ mb: 1 }}>
                    {HELP_DIALOG_MESSAGE[0]}{" "}
                    <Link
                        href={`mailto:${MDV_EMAIL}?subject=MDV%20Feedback%20or%20Issue`}
                        sx={{
                            textDecoration: "none",
                            "&.MuiLink-root": { color: "info.main" },
                        }}
                    >
                        {MDV_EMAIL}
                    </Link>
                    .
                </Typography>
                <Typography  variant="body1" sx={{mt: 1}}>
                    {HELP_DIALOG_MESSAGE[1]}
                </Typography>
                <Typography sx={{mt: 3, mb: 1 }} variant="body1">
                    {HELP_DIALOG_MESSAGE[2]}{" "}
                    <Link
                        href={MDV_DOC_WEBSITE}
                        target="_blank"
                        rel="noopener noreferrer"
                        sx={{
                            textDecoration: "none",
                            "&.MuiLink-root": { color: "info.main" },
                        }}
                    >
                        {MDV_DOC_WEBSITE}
                    </Link>
                </Typography>
            </DialogContent>
            <DialogActions>
                <Button onClick={onClose} color="error">
                    Close
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default HelpDialog;
