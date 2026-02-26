import {
    Alert,
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    TextField,
} from "@mui/material";
import { useEffect, useState } from "react";
import { DEFAULT_GATE_COLOR } from "../gates/gateUtils";
import { hexToRgb, rgbToHex } from "@/utilities/Utilities";

export type GateNameDialogType = {
    open: boolean;
    onClose: () => void;
    onSaveGate: (gateName: string, color?: [number, number, number]) => void;
    name?: string;
    isEditName?: boolean;
};

const GateNameDialog = ({ open, onClose, onSaveGate, name, isEditName }: GateNameDialogType) => {
    const [gateName, setGateName] = useState(name ? name : "");
    const [color, setColor] = useState<[number, number, number]>(DEFAULT_GATE_COLOR);
    const [error, setError] = useState("");

    useEffect(() => {
        if (open) {
            setGateName(name ?? "");
            setColor(DEFAULT_GATE_COLOR);
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
            if (isEditName) {
                onSaveGate(gateName);
            } else {
                onSaveGate(gateName, color);
            }
            handleClose();
        } catch (err) {
            const errorMsg = err instanceof Error ? err.message : "Something went wrong while saving the gate";
            setError(errorMsg);
        }
    };

    const handleClose = () => {
        setGateName("");
        setColor(DEFAULT_GATE_COLOR);
        setError("");
        onClose();
    };

    return (
        <Dialog open={open} onClose={handleClose} maxWidth="sm" fullWidth>
            <DialogTitle>Gate Name</DialogTitle>
            <DialogContent dividers>
                <Box sx={{ display: "flex", alignItems: "center", gap: 2, p: 2 }}>
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
                    />
                    {!isEditName ? (
                        <input
                            type="color"
                            value={rgbToHex(color)}
                            onChange={(e) => setColor(hexToRgb(e.target.value))}
                            style={{
                                width: 35,
                                height: 35,
                                borderRadius: 4,
                                cursor: "pointer",
                            }}
                            aria-label="Gate color"
                        />
                        ): null
                    }  
                </Box>

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
