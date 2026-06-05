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

type RenameTableColumnDialogProps = {
    open: boolean;
    columnField: string | null;
    initialName: string;
    onClose: () => void;
    onSubmit: (args: { columnField: string; newName: string }) => void;
};

const RenameTableColumnDialog = ({
    open,
    columnField,
    initialName,
    onClose,
    onSubmit,
}: RenameTableColumnDialogProps) => {
    const [name, setName] = useState(initialName);
    const [nameError, setNameError] = useState<string | null>(null);

    useEffect(() => {
        if (!open) {
            return;
        }
        setName(initialName);
        setNameError(null);
    }, [open, initialName]);

    const submit = () => {
        if (!columnField) {
            return;
        }
        const trimmedName = name.trim();
        if (!trimmedName) {
            setNameError("Column name is required");
            return;
        }
        if (trimmedName.toLocaleLowerCase() === initialName.trim().toLocaleLowerCase()) {
            onClose();
            return;
        }
        onSubmit({
            columnField,
            newName: trimmedName,
        });
    };

    return (
        <Dialog open={open} onClose={onClose} fullWidth maxWidth="sm">
            <DialogTitle>
                {initialName ? `Rename "${initialName}"` : "Rename Column"}
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers>
                <Stack spacing={2} sx={{ pt: 0.5 }}>
                    <Typography variant="subtitle2" color="text.secondary">
                        Change the display name for this column. The underlying field id will stay the same.
                    </Typography>
                    <TextField
                        autoFocus
                        label="Column Name"
                        value={name}
                        onChange={(event) => {
                            setName(event.target.value);
                            setNameError(null);
                        }}
                        onKeyDown={(event) => {
                            if (event.key === "Enter") {
                                event.preventDefault();
                                submit();
                            }
                        }}
                        error={!!nameError}
                        helperText={nameError}
                        fullWidth
                    />
                </Stack>
            </DialogContent>
            <DialogActions sx={{ flexWrap: "wrap", gap: 1 }}>
                <Button variant="contained" onClick={submit}>
                    Rename
                </Button>
                <Button color="error" onClick={onClose}>
                    Cancel
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default RenameTableColumnDialog;
