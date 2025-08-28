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
    TextField,
    Typography,
} from "@mui/material";
import { useState } from "react";
import type React from "react";

export interface AnndataConflictDialogProps {
    open: boolean;
    onClose: () => void;
    onConfirm: (prefix: string) => void;
    onError?: (error: { message: string; status?: number }) => void;
}

const AnndataConflictDialog: React.FC<AnndataConflictDialogProps> = ({
    open,
    onClose,
    onConfirm,
    onError,
}) => {
    const [prefix, setPrefix] = useState("");
    const [error, setError] = useState("");

    const isValidPrefix =
        prefix.trim() && /^[A-Za-z][A-Za-z0-9_]*$/.test(prefix);

    const handleConfirm = () => {
        if (!prefix.trim()) {
            const errorMessage = "Prefix cannot be empty";
            setError(errorMessage);
            onError?.({ message: errorMessage, status: 400 });
            return;
        }
        if (!/^[A-Za-z][A-Za-z0-9_]*$/.test(prefix)) {
            const errorMessage =
                "Prefix must start with a letter and contain only letters, numbers, and underscores";
            setError(errorMessage);
            onError?.({ message: errorMessage, status: 400 });
            return;
        }
        onConfirm(prefix);
    };

    const handlePrefixChange = (event: React.ChangeEvent<HTMLInputElement>) => {
        setPrefix(event.target.value);
        setError("");
    };

    return (
        <Dialog
            open={open}
            onClose={onClose}
            PaperComponent={Paper}
            hideBackdrop={true}
            PaperProps={{
                elevation: 24,
                className: "rounded-lg",
            }}
            sx={{
                pointerEvents: "none",
                "& .MuiPaper-root": {
                    pointerEvents: "auto",
                },
            }}
        >
            <DialogTitle className="flex items-center gap-3 bg-blue-50 dark:bg-blue-900 border-b px-6 py-4">
                <MergeTypeIcon className="text-blue-600 dark:text-blue-400" />
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
                        An AnnData file already exists in this project. To
                        combine the files, please specify a prefix for the new
                        datasources to prevent naming conflicts:
                    </Typography>

                    <TextField
                        fullWidth
                        label="Prefix for new datasources"
                        value={prefix}
                        onChange={handlePrefixChange}
                        error={!!error}
                        helperText={
                            error ||
                            "This prefix will be added to all new datasource names"
                        }
                        placeholder="Enter prefix"
                        className="mt-4"
                        autoFocus
                    />

                    <Box className="bg-gray-50 dark:bg-gray-900 p-4 rounded-md font-mono text-sm">
                        <div>existing_cells → cells</div>
                        <div>new_cells → {prefix || "(prefix)"}_cells</div>
                    </Box>

                    <Alert severity="info" className="mt-4">
                        This will preserve your existing data, while adding new
                        datasources with your chosen prefix (columns will not be
                        merged, these will be distinct, separate, datasources).
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
                    onClick={handleConfirm}
                    variant="contained"
                    color="primary"
                    className="bg-blue-600 hover:bg-blue-700"
                    disabled={!isValidPrefix}
                >
                    Combine Sources
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default AnndataConflictDialog;
