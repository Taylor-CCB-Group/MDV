import { Alert, Button, Dialog, DialogActions, DialogContent, DialogTitle, TextField } from "@mui/material";
import { useEffect, useState } from "react";

export type GateNameDialogType = {
    open: boolean;
    onClose: () => void;
    onSaveGate: (gateName: string) => void;
    name?: string;
};

const GateNameDialog = ({ open, onClose, onSaveGate, name }: GateNameDialogType) => {
    const [gateName, setGateName] = useState(name ? name : "");
    const [error, setError] = useState("");

    useEffect(() => {
        if (open) {
            setGateName(name ?? "");
            setError("");
        }
    }, [open, name]);

    const handleSave = () => {
        // Validate name
        if (!gateName.trim()) {
            setError("Gate name cannot be empty");
            return;
        }

        try {
            onSaveGate(gateName);
            handleClose();
        } catch (err) {
            const errorMsg = err instanceof Error ? err.message : "Something went wrong while saving the gate";
            setError(errorMsg);
        }
    };

    const handleClose = () => {
        setGateName("");
        setError("");
        onClose();
    };

    return (
        <Dialog open={open} onClose={handleClose} maxWidth="sm" fullWidth>
            <DialogTitle>Gate Name</DialogTitle>
            <DialogContent dividers>
                <TextField
                    autoFocus
                    fullWidth
                    label="Gate name"
                    placeholder="e.g., CD4+ T cells"
                    value={gateName}
                    onChange={(e) => {
                        setGateName(e.target.value);
                        setError("");
                    }}
                    onKeyDown={(e) => e.key === "Enter" && handleSave()}
                    sx={{ mt: 1 }}
                />

                {error && (
                    <Alert
                        severity="error"
                        sx={{
                            mt: 2,
                            backgroundColor: "var(--background_color_error)",
                        }}
                    >
                        {error}
                    </Alert>
                )}
            </DialogContent>
            <DialogActions>
                <Button onClick={handleClose}>Cancel</Button>
                <Button onClick={handleSave} variant="contained">
                    Save Gate
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default GateNameDialog;
