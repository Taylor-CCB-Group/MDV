import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";
import {
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Stack,
    TextField,
    Typography,
} from "@mui/material";
import { useEffect, useState } from "react";

export type BulkEditAction =
    | "fill-all"
    | "fill-empty";

type BulkEditColumnDialogProps = {
    open: boolean;
    columnName: string | null;
    onClose: () => void;
    onSubmit: (args: {
        action: BulkEditAction;
        columnName: string;
        value: string;
    }) => void;
};

const BulkEditColumnDialog = ({
    open,
    columnName,
    onClose,
    onSubmit,
}: BulkEditColumnDialogProps) => {
    const [value, setValue] = useState("");

    useEffect(() => {
        if (!open) {
            return;
        }
        setValue("");
    }, [open, columnName]);

    const submit = (action: BulkEditAction) => {
        if (!columnName) {
            return;
        }
        onSubmit({
            action,
            columnName,
            value,
        });
    };

    return (
        <Dialog open={open} onClose={onClose} fullWidth maxWidth="sm">
            <DialogTitle>
                {columnName ? `Bulk Edit in "${columnName}"` : "Bulk Edit Column"}
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers>
                <Stack spacing={2} sx={{ pt: 0.5 }}>
                    <Typography variant="subtitle2" color="text.secondary">
                        Use this dialog to fill all visible cells or only the empty cells in the selected column.
                    </Typography>
                    <TextField
                        label="Value"
                        value={value}
                        onChange={(event) => setValue(event.target.value)}
                        fullWidth
                    />
                </Stack>
            </DialogContent>
            <DialogActions sx={{ flexWrap: "wrap", gap: 1 }}>
                <Button variant="contained" onClick={() => submit("fill-all")}>
                    Fill All Cells
                </Button>
                <Button variant="contained" onClick={() => submit("fill-empty")}>
                    Fill Empty Cells
                </Button>
                <Button color="error" onClick={onClose}>
                    Cancel
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default BulkEditColumnDialog;
