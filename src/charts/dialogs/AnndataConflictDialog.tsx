import MergeTypeIcon from "@mui/icons-material/MergeType";
import {
    Alert,
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Paper,
    Typography,
} from "@mui/material";
import type React from "react";

interface AnndataConflictDialogProps {
    open: boolean;
    onClose: () => void;
    onConfirm: () => void;
}

const AnndataConflictDialog: React.FC<AnndataConflictDialogProps> = ({
    open,
    onClose,
    onConfirm,
}) => {
    return (
        <Dialog
            open={open}
            onClose={onClose}
            PaperComponent={Paper}
            maxWidth="sm"
            fullWidth
            PaperProps={{
                elevation: 24,
                className: "rounded-lg",
            }}
        >
            <DialogTitle className="flex items-center gap-3 bg-blue-50 dark:bg-blue-900 border-b px-6 py-4">
                <Typography
                    variant="h6"
                    component="span"
                    className="font-semibold"
                >
                    Existing AnnData File Found
                </Typography>
            </DialogTitle>

            <DialogContent className="p-6">
                <Box className="space-y-4 mt-4">
                    <Typography variant="body1">
                        An AnnData file, along with its datasources, already
                        exists in this project. If you choose to add a new
                        AnnData file, all datasources will automatically be
                        assigned unique prefixes to prevent naming conflicts:
                    </Typography>

                    <Box className="bg-gray-50 dark:bg-gray-900 p-4 rounded-md font-mono text-sm">
                        <div>existing_cells → A_cells</div>
                        <div>new_cells → B_cells</div>
                    </Box>

                    <Alert severity="info" className="mt-4">
                        This will preserve your existing data while adding the
                        new sources.
                    </Alert>
                </Box>
            </DialogContent>

            <DialogActions className="bg-gray-50 dark:bg-gray-900 border-t p-4 gap-2">
                <Button
                    onClick={onClose}
                    variant="outlined"
                    color="inherit"
                    className="hover:bg-gray-100 dark:hover:bg-gray-800"
                >
                    Cancel
                </Button>
                <Button
                    onClick={onConfirm}
                    variant="contained"
                    color="primary"
                    className="bg-blue-600 hover:bg-blue-700"
                    autoFocus
                >
                    Combine Sources
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default AnndataConflictDialog;
