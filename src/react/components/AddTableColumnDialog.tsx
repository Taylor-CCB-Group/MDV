import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";
import type { DataType } from "@/charts/charts";
import {
    Autocomplete,
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

type CloneableColumn = {
    field: string;
    name: string;
    datatype: DataType;
    stringLength?: number;
    delimiter?: string;
};

export type AddColumnParams = {
    name: string;
    datatype: DataType;
    cloneColumn?: string | null;
    sourceColumns?: string[] | null;
    mode?: "empty" | "clone" | "compound";
    position?: number | null;
    stringLength?: number | null;
    delimiter?: string | null;
};

type AddColumnMode = "empty" | "clone" | "compound";

type AddTableColumnDialogProps = {
    open: boolean;
    cloneableColumns: CloneableColumn[];
    defaultPosition: number;
    onClose: () => void;
    onSubmit: (args: AddColumnParams) => void;
};

//! We currently don't support int32 datatype column creation
//! int32 cannot truly represent NaN, so empty int32 cells currently coerce during storage and need a real missing-value strategy later.
const DATATYPE_OPTIONS: DataType[] = [
    "text",
    "text16",
    "integer",
    "double",
    "unique",
    "multitext",
];
const COMPOUND_DATATYPE_OPTIONS: Array<"text" | "text16"> = ["text", "text16"];
const ADD_COLUMN_MODE_OPTIONS: AddColumnMode[] = ["empty", "clone", "compound"];

const AddTableColumnDialog = ({
    open,
    cloneableColumns,
    defaultPosition,
    onClose,
    onSubmit,
}: AddTableColumnDialogProps) => {
    const [name, setName] = useState("");
    const [mode, setMode] = useState<AddColumnMode>("empty");
    const [datatype, setDatatype] = useState<DataType>("text");
    const [cloneColumn, setCloneColumn] = useState("");
    const [sourceColumns, setSourceColumns] = useState<string[]>([]);
    const [position, setPosition] = useState(String(defaultPosition));
    const [stringLength, setStringLength] = useState("");
    const [delimiter, setDelimiter] = useState(",");
    const [nameError, setNameError] = useState<string | null>(null);
    const [sourceColumnError, setSourceColumnError] = useState<string | null>(null);

    const selectedCloneColumn = cloneableColumns.find((column) => column.field === cloneColumn);
    const compoundSourceColumns = cloneableColumns.filter(
        (column) => column.datatype === "text" || column.datatype === "text16",
    );
    const selectedSourceColumns = compoundSourceColumns.filter((column) =>
        sourceColumns.includes(column.field),
    );
    const selectedDatatype =
        mode === "clone"
            ? selectedCloneColumn?.datatype ?? "text"
            : mode === "compound"
              ? (datatype === "text16" ? "text16" : "text")
              : datatype;

    // Reset values when the dialog is opened to avoid stale values
    useEffect(() => {
        if (!open) {
            return;
        }
        setName("");
        setMode("empty");
        setDatatype("text");
        setCloneColumn(cloneableColumns[0]?.field ?? "");
        setSourceColumns([]);
        setPosition(String(defaultPosition));
        setStringLength("");
        setDelimiter("_");
        setNameError(null);
        setSourceColumnError(null);
    }, [open, cloneableColumns, defaultPosition]);

    const handleSubmit = () => {
        const trimmedName = name.trim();
        if (!trimmedName) {
                setNameError("Column name is required");
            return;
        }
        if (mode === "compound" && sourceColumns.length < 2) {
            setSourceColumnError("Select at least 2 source columns");
            return;
        }
        const parsedPosition = Number.parseInt(position, 10);
        const parsedStringLength = Number.parseInt(stringLength, 10);
        onSubmit({
            name: trimmedName,
            mode,
            datatype: selectedDatatype,
            cloneColumn: mode === "clone" && cloneColumn ? cloneColumn : null,
            sourceColumns: mode === "compound" ? sourceColumns : null,
            position: Number.isNaN(parsedPosition) ? null : parsedPosition,
            stringLength:
                selectedDatatype === "unique" || selectedDatatype === "multitext"
                    ? Number.isNaN(parsedStringLength)
                        ? null
                        : parsedStringLength
                    : null,
            delimiter:
                selectedDatatype === "multitext" || mode === "compound"
                    ? delimiter
                    : null,
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
                        onChange={(event) => {
                                setName(event.target.value);
                                setNameError(null);
                            }}
                        onKeyDown={(event) => {
                            if (event.key === "Enter") {
                                event.preventDefault();
                                handleSubmit();
                            }
                        }}
                            error={!!nameError}
                            helperText={nameError}
                        fullWidth
                    />

                    <Autocomplete
                        options={ADD_COLUMN_MODE_OPTIONS}
                        value={mode}
                        onChange={(_, value) => setMode(value ?? "empty")}
                        disableClearable
                        getOptionLabel={(option) =>
                            option === "empty"
                                ? "Empty"
                                : option === "clone"
                                  ? "Clone"
                                  : "Compound"
                        }
                        renderInput={(params) => (
                            <TextField {...params} label="Mode" fullWidth />
                        )}
                    />

                    {mode === "clone" ? (
                        <>
                            <Autocomplete
                                options={cloneableColumns}
                                value={selectedCloneColumn}
                                onChange={(_, value) => {
                                    if (value) {
                                        setCloneColumn(value.field);
                                    }
                                }}
                                isOptionEqualToValue={(option, value) => option.field === value.field}
                                disableClearable
                                getOptionLabel={(option) =>
                                    typeof option === "string"
                                        ? option
                                        : option.name
                                }
                                renderOption={(props, option) => (
                                    <li {...props}>
                                            <Typography variant="body1">{option.name}</Typography>
                                            <Typography
                                                variant="body2"
                                                color="textDisabled"
                                                sx={{
                                                    ml: 1,
                                                    fontStyle: "italic",
                                                }}
                                            >
                                                ({option.datatype})
                                            </Typography>
                                    </li>
                                )}
                                disabled={cloneableColumns.length === 0}
                                renderInput={(params) => (
                                    <TextField
                                        {...params}
                                        label="Column to clone"
                                        placeholder="Search columns"
                                        fullWidth
                                    />
                                )}
                            />

                            <TextField
                                label="Inherited Datatype"
                                value={selectedDatatype}
                                fullWidth
                                disabled
                                slotProps={{
                                    input: {
                                        readOnly: true,
                                    }
                                }}
                            />
                        </>
                    ) : mode === "compound" ? (
                        <>
                            <Autocomplete
                                options={compoundSourceColumns}
                                multiple
                                value={selectedSourceColumns}
                                onChange={(_, values) => {
                                    setSourceColumns(values.map((value) => value.field));
                                    setSourceColumnError(null);
                                }}
                                getOptionLabel={(option) => option.name}
                                isOptionEqualToValue={(option, value) => option.field === value.field}
                                renderInput={(params) => (
                                    <TextField
                                        {...params}
                                        label="Source Columns"
                                        error={!!sourceColumnError}
                                        helperText={sourceColumnError}
                                        fullWidth
                                    />
                                )}
                            />
                            <Autocomplete
                                options={COMPOUND_DATATYPE_OPTIONS}
                                value={selectedDatatype}
                                onChange={(_, value) => setDatatype(value ?? "text")}
                                getOptionLabel={(option) => option}
                                renderInput={(params) => (
                                    <TextField {...params} label="Datatype" fullWidth />
                                )}
                            />
                            <TextField
                                label="Separator"
                                value={delimiter}
                                onChange={(event) => setDelimiter(event.target.value)}
                                fullWidth
                            />
                        </>
                    ) : (
                        <Autocomplete
                            options={DATATYPE_OPTIONS}
                            value={datatype}
                            onChange={(_, value) => setDatatype(value ?? "text")}
                            getOptionLabel={(option) => option}
                            renderInput={(params) => (
                                <TextField
                                    {...params}
                                    label="Datatype"
                                    placeholder="Search datatype"
                                    fullWidth
                                />
                            )}
                        />
                    )}

                    {selectedDatatype === "unique" && mode === "empty" && (
                        <TextField
                            label="Max Length"
                            value={stringLength}
                            onChange={(event) => setStringLength(event.target.value)}
                            type="number"
                            slotProps={{
                                htmlInput: {
                                    min: 1,
                                }
                            }}
                            fullWidth
                        />
                    )}

                    {selectedDatatype === "multitext" && mode === "empty" && (
                        <>
                            <TextField
                                label="Capacity"
                                value={stringLength}
                                onChange={(event) => setStringLength(event.target.value)}
                                type="number"
                                slotProps={{
                                    htmlInput: {
                                        min: 1,
                                    }
                                }}
                                fullWidth
                            />
                            <TextField
                                label="Delimiter"
                                value={delimiter}
                                onChange={(event) => setDelimiter(event.target.value)}
                                fullWidth
                            />
                        </>
                    )}

                    <TextField
                        label="Position"
                        value={position}
                        onChange={(event) => setPosition(event.target.value)}
                        type="number"
                        slotProps={{
                            htmlInput: {
                                min: 1,
                            }
                        }}
                        fullWidth
                    />

                    <Typography variant="body2" color="text.secondary">
                        New columns are editable. Compound source columns are limited to text/text16, and a requested text output auto-upgrades to text16 if distinct values exceed 256.
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
