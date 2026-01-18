import { Alert, Button, Dialog, DialogActions, DialogContent, DialogTitle, TextField, Typography } from "@mui/material";
import { useState } from "react";
import { observer } from 'mobx-react-lite';

export type GateNameDialogType = {
    open: boolean;
    onClose: () => void;
    onSaveGate: (gateName: string) => void;
};

const GateNameDialog = observer(function GateNameDialog({ 
    open, 
    onClose,
    onSaveGate,
}: GateNameDialogType) {
    const [gateName, setGateName] = useState('');
    const [error, setError] = useState('');

    const handleSave = () => {
        // Validate name
        if (!gateName.trim()) {
            setError('Gate name cannot be empty');
            return;
        }

        try {
            onSaveGate(gateName);
        } catch (err) {
            const errorMsg = err instanceof Error ? err.message : "Something went wrong while saving the gate";
            setError(errorMsg);
        }
        
        // Close
        handleClose();
    };

    const handleClose = () => {
        setGateName("");
        setError("");
        onClose();
    };

    return (
        <Dialog open={open} onClose={handleClose} maxWidth="sm" fullWidth>
            <DialogTitle>Create Gate</DialogTitle>
            <DialogContent>
                <Typography variant="caption" color="textSecondary" paragraph>
                    Gates define cell populations. This gate will appear on all charts with matching axes.
                </Typography>
                
                {error && <Alert severity="error" sx={{ mb: 2 }}>{error}</Alert>}
                
                <TextField
                    autoFocus
                    fullWidth
                    label="Gate name"
                    placeholder="e.g., CD4+ T cells"
                    value={gateName}
                    onChange={(e) => {
                        setGateName(e.target.value);
                        setError('');
                    }}
                    onKeyPress={(e) => e.key === 'Enter' && handleSave()}
                    sx={{ mt: 1 }}
                />
            </DialogContent>
            <DialogActions>
                <Button onClick={handleClose}>Cancel</Button>
                <Button onClick={handleSave} variant="contained">Save Gate</Button>
            </DialogActions>
        </Dialog>
    );
});

export default GateNameDialog;

