import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";
import {
    Button,
    Checkbox,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    FormControlLabel,
    MenuItem,
    Select,
    Stack,
    TextField,
    Typography,
} from "@mui/material";
import { useEffect, useState } from "react";

type CloneableColumn = {
    field: string;
    name: string;
};

type AddTableColumnDialogProps = {
    open: boolean;
    cloneableColumns: CloneableColumn[];
    defaultPosition: number;
    onClose: () => void;
    onSubmit: (args: {
        name: string;
        cloneColumn: string | null;
        position: number | null;
    }) => void;
};

const AddTableColumnDialog = ({
    open,
    cloneableColumns,
    defaultPosition,
    onClose,
    onSubmit,
}: AddTableColumnDialogProps) => {
    const [name, setName] = useState("");
    const [copyExistingColumn, setCopyExistingColumn] = useState(false);
    const [cloneColumn, setCloneColumn] = useState("");
    const [position, setPosition] = useState(String(defaultPosition));

    useEffect(() => {
        if (!open) {
            return;
        }
        setName("");
        setCopyExistingColumn(false);
        setCloneColumn(cloneableColumns[0]?.field ?? "");
        setPosition(String(defaultPosition));
    }, [open, cloneableColumns, defaultPosition]);

    const handleSubmit = () => {
        const parsedPosition = Number.parseInt(position, 10);
        onSubmit({
            name,
            cloneColumn: copyExistingColumn && cloneColumn ? cloneColumn : null,
            position: Number.isNaN(parsedPosition) ? null : parsedPosition,
        });
    };

    return (
        <Dialog open={open} onClose={onClose} fullWidth maxWidth="xs">
            <DialogTitle>
                Add Column
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers>
                <Stack spacing={2} sx={{ pt: 0.5 }}>
                    <TextField
                        autoFocus
                        label="Column Name"
                        value={name}
                        onChange={(event) => setName(event.target.value)}
                        onKeyDown={(event) => {
                            if (event.key === "Enter") {
                                event.preventDefault();
                                handleSubmit();
                            }
                        }}
                        fullWidth
                    />

                    <FormControlLabel
                        control={
                            <Checkbox
                                checked={copyExistingColumn}
                                onChange={(event) => setCopyExistingColumn(event.target.checked)}
                                disabled={cloneableColumns.length === 0}
                            />
                        }
                        label="Copy Existing Column"
                    />

                    {copyExistingColumn && (
                        <Select
                            value={cloneColumn}
                            onChange={(event) => setCloneColumn(String(event.target.value))}
                            disabled={!copyExistingColumn || cloneableColumns.length === 0}
                            fullWidth
                            displayEmpty
                        >
                            {cloneableColumns.length === 0 ? (
                                <MenuItem value="" disabled>
                                    No text columns available
                                </MenuItem>
                            ) : (
                                cloneableColumns.map((column) => (
                                    <MenuItem key={column.field} value={column.field}>
                                        {column.name}
                                    </MenuItem>
                                ))
                            )}
                        </Select>
                    )}

                    <TextField
                        label="Position"
                        value={position}
                        onChange={(event) => setPosition(event.target.value)}
                        type="number"
                        inputProps={{ min: 1 }}
                        fullWidth
                    />

                    <Typography variant="body2" color="text.secondary">
                        New columns are editable and can be used for annotations.
                    </Typography>
                </Stack>
            </DialogContent>
            <DialogActions>
                <Button onClick={handleSubmit} variant="contained">
                    Add
                </Button>
                <Button onClick={onClose} color="error">Cancel</Button>
            </DialogActions>
        </Dialog>
    );
};

export default AddTableColumnDialog;
