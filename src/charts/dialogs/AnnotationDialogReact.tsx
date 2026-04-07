import { type PropsWithChildren, useEffect, useMemo, useState } from "react";
import type DataStore from "../../datastore/DataStore.js";
import TagModel, { type TagSelectionScope } from "../../table/TagModel";
import { BaseDialog } from "../../utilities/Dialog.js";
import { createMdvPortal } from "@/react/react_utils";
import { observer } from "mobx-react-lite";
import {
    Autocomplete,
    Box,
    Button,
    Checkbox,
    Chip,
    Divider,
    FormControl,
    FormControlLabel,
    Paper,
    Radio,
    RadioGroup,
    Stack,
    TextField,
    Typography,
} from "@mui/material";
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";

const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;

function TagView({
    dataStore,
    columnName,
    selectionScope,
}: {
    dataStore: DataStore;
    columnName: string;
    selectionScope: TagSelectionScope;
}) {
    const [tagModel, setTagModel] = useState<TagModel | null>(null);
    const [loadError, setLoadError] = useState<string | null>(null);
    const [tagList, setTagList] = useState<Set<string>>(new Set());
    const [tagsInSelection, setTagsInSelection] = useState<Set<string>>(new Set());
    const [draftTag, setDraftTag] = useState("");
    const [selectionCount, setSelectionCount] = useState(0);

    useEffect(() => {
        let cancelled = false;

        const syncFromModel = (model: TagModel) => {
            setTagList(model.getTags());
            setTagsInSelection(model.getTagsInSelection());
            setSelectionCount(model.getSelectionLength());
        };

        setTagModel(null);
        setLoadError(null);
        TagModel.create(dataStore, columnName, selectionScope)
            .then((model) => {
                if (cancelled) {
                    return;
                }
                setTagModel(model);
                syncFromModel(model);
            })
            .catch((error) => {
                if (cancelled) {
                    return;
                }
                const message =
                    error instanceof Error ? error.message : String(error);
                setLoadError(message);
            });

        return () => {
            cancelled = true;
        };
    }, [columnName, dataStore, selectionScope]);

    useEffect(() => {
        if (!tagModel) {
            return;
        }

        const syncFromModel = () => {
            tagModel.setSelectionScope(selectionScope);
            setTagList(tagModel.getTags());
            setTagsInSelection(tagModel.getTagsInSelection());
            setSelectionCount(tagModel.getSelectionLength());
        };

        const listener = tagModel.addListener(syncFromModel);
        syncFromModel();
        return () => {
            tagModel.removeListener(listener);
        };
    }, [selectionScope, tagModel]);

    const availableTags = useMemo(
        () =>
            Array.from(tagList).sort((left, right) => left.localeCompare(right)),
        [tagList],
    );
    const selectedTags = useMemo(
        () =>
            Array.from(tagsInSelection).sort((left, right) =>
                left.localeCompare(right),
            ),
        [tagsInSelection],
    );

    if (loadError) {
        return (
            <Typography color="error" variant="body2">
                Could not load annotation column: {loadError}
            </Typography>
        );
    }

    if (!tagModel) {
        return <Typography variant="body2">Loading annotations...</Typography>;
    }

    const addTag = (rawTag: string) => {
        const nextTag = rawTag.trim();
        if (!nextTag) {
            return;
        }
        tagModel.setTag(nextTag, true);
        setDraftTag("");
    };

    return (
        <Paper
            sx={{
                border: "1px solid",
                borderColor: "divider",
                p: 2,
            }}
            variant="outlined"
        >
            <Stack spacing={2.5}>
                <Box>
                <Typography variant="h6">
                    Annotating {selectionCount} {selectionScope} row{selectionCount === 1 ? "" : "s"}
                </Typography>
                <Typography color="text.secondary" variant="body2">
                    Column: {columnName}
                </Typography>
                </Box>

                <Stack direction="row" spacing={1}>
                <Autocomplete
                    freeSolo
                    fullWidth
                    size="small"
                    options={availableTags}
                    inputValue={draftTag}
                    onInputChange={(_, value) => {
                        setDraftTag(value);
                    }}
                    onChange={(_, value) => {
                        if (typeof value === "string") {
                            addTag(value);
                        }
                    }}
                    renderInput={(params) => {
                        const { key, ...fieldProps } = params as typeof params & {
                            key: string;
                        };
                        return (
                            <TextField
                                key={key}
                                {...fieldProps}
                                label="Add annotation item"
                                placeholder="Type a value and press Enter"
                                onKeyDown={(event) => {
                                    if (event.key !== "Enter") {
                                        return;
                                    }
                                    event.preventDefault();
                                    addTag(draftTag);
                                }}
                            />
                        );
                    }}
                />
                <Button
                    onClick={() => addTag(draftTag)}
                    sx={{ minWidth: 96 }}
                    variant="contained"
                >
                    Add
                </Button>
                </Stack>

                <Box>
                <Typography gutterBottom variant="subtitle2">
                    Present on selected rows
                </Typography>
                {selectedTags.length === 0 ? (
                    <Typography color="text.secondary" variant="body2">
                        No annotation items are present on the current selection.
                    </Typography>
                ) : (
                    <Box sx={{ display: "flex", flexWrap: "wrap", gap: 1 }}>
                        {selectedTags.map((tag) => {
                            const entireSelection = tagModel.entireSelectionHasTag(tag);
                            return (
                                <Chip
                                    color={entireSelection ? "primary" : "default"}
                                    key={tag}
                                    label={tag}
                                    onDelete={() => tagModel.setTag(tag, false)}
                                    variant={entireSelection ? "filled" : "outlined"}
                                />
                            );
                        })}
                    </Box>
                )}
                </Box>

                <Divider />

                <Box>
                <Typography gutterBottom variant="subtitle2">
                    Available items
                </Typography>
                {availableTags.length === 0 ? (
                    <Typography color="text.secondary" variant="body2">
                        Add the first annotation item for this column above.
                    </Typography>
                ) : (
                    <Paper
                        sx={{
                            maxHeight: 240,
                            overflowY: "auto",
                            border: "1px solid",
                            borderColor: "divider",
                        }}
                        variant="outlined"
                    >
                        {availableTags.map((tag) => {
                            const entireSelection = tagModel.entireSelectionHasTag(tag);
                            const inSelection = tagsInSelection.has(tag);
                            const indeterminate = inSelection && !entireSelection;
                            return (
                                <Box
                                    key={tag}
                                    onClick={() =>
                                        tagModel.setTag(tag, indeterminate ? true : !entireSelection)
                                    }
                                    sx={{
                                        alignItems: "center",
                                        cursor: "pointer",
                                        display: "flex",
                                        gap: 1,
                                        px: 1.5,
                                        py: 1,
                                        "&:hover": {
                                            backgroundColor: "action.hover",
                                        },
                                    }}
                                >
                                    <Checkbox
                                        checked={entireSelection}
                                        checkedIcon={checkedIcon}
                                        icon={icon}
                                        indeterminate={indeterminate}
                                        onChange={() => undefined}
                                        sx={{ p: 0.5 }}
                                    />
                                    <Box sx={{ minWidth: 0 }}>
                                        <Typography noWrap variant="body2">
                                            {tag}
                                        </Typography>
                                        <Typography color="text.secondary" variant="caption">
                                            {entireSelection
                                                ? "Present on all selected rows"
                                                : inSelection
                                                  ? "Present on some selected rows"
                                                  : "Not present on the current selection"}
                                        </Typography>
                                    </Box>
                                </Box>
                            );
                        })}
                    </Paper>
                )}
                </Box>
            </Stack>
        </Paper>
    );
}

function SectionCard({
    title,
    description,
    children,
}: PropsWithChildren<{
    title: string;
    description: string;
}>) {
    return (
        <Paper
            sx={{
                border: "1px solid",
                borderColor: "divider",
                p: 2,
            }}
            variant="outlined"
        >
            <Stack spacing={2}>
                <Box>
                    <Typography variant="h6">{title}</Typography>
                    <Typography color="text.secondary" variant="body2">
                        {description}
                    </Typography>
                </Box>
                {children}
            </Stack>
        </Paper>
    );
}

const AnnotationDialogComponent = observer(
    ({ dataStore }: { dataStore: DataStore }) => {
        const columns = useMemo(
            () =>
                dataStore.columns
                    .filter((column) => column.datatype === "multitext")
                    .map((column) => column.field || column.name)
                    .sort((left, right) => left.localeCompare(right)),
            [dataStore],
        );
        const [selectedColumn, setSelectedColumn] = useState(
            () => dataStore.columnIndex.__tags ? "__tags" : columns[0] || "",
        );
        const [columnInput, setColumnInput] = useState(selectedColumn);
        const [selectionScope, setSelectionScope] =
            useState<TagSelectionScope>("filtered");

        const getNameState = (rawValue: string) => {
            const nextValue = rawValue.trim();
            if (!nextValue) {
                return "empty" as const;
            }
            const existingColumn = dataStore.columnIndex[nextValue];
            return existingColumn && existingColumn.datatype !== "multitext"
                ? ("clash" as const)
                : ("ok" as const);
        };

        const nameState = getNameState(columnInput);

        const applyColumnChoice = (rawValue: string) => {
            const nextValue = rawValue.trim();
            if (getNameState(nextValue) !== "ok") {
                return;
            }
            setSelectedColumn(nextValue);
            setColumnInput(nextValue);
        };

        return (
            <Stack spacing={2}>
                <SectionCard
                    title="Annotation Setup"
                    description="Choose the multitext column to edit, then pick whether tag actions should apply to filtered rows or the current highlight set."
                >
                    <Stack spacing={2}>
                        <Stack direction="row" spacing={1}>
                            <Autocomplete
                                freeSolo
                                fullWidth
                                size="small"
                                inputValue={columnInput}
                                onChange={(_, value) => {
                                    if (typeof value === "string") {
                                        setColumnInput(value);
                                        applyColumnChoice(value);
                                    }
                                }}
                                onInputChange={(_, value) => {
                                    setColumnInput(value);
                                }}
                                options={columns}
                                renderInput={(params) => {
                                    const { key, ...fieldProps } = params as typeof params & {
                                        key: string;
                                    };
                                    return (
                                        <TextField
                                            key={key}
                                            {...fieldProps}
                                            error={nameState === "clash"}
                                            helperText={
                                                nameState === "clash"
                                                    ? "A non-multitext column already uses that name."
                                                    : "Names become new multitext annotation columns when needed."
                                            }
                                            label="Column name"
                                            onKeyDown={(event) => {
                                                if (event.key !== "Enter") {
                                                    return;
                                                }
                                                event.preventDefault();
                                                applyColumnChoice(columnInput);
                                            }}
                                        />
                                    );
                                }}
                            />
                            <Button
                                disabled={nameState !== "ok"}
                                onClick={() => applyColumnChoice(columnInput)}
                                sx={{ minWidth: 112 }}
                                variant="contained"
                            >
                                Use Column
                            </Button>
                        </Stack>

                        <FormControl>
                            <Typography gutterBottom variant="subtitle2">
                                Apply actions to
                            </Typography>
                            <RadioGroup
                                row
                                onChange={(_, value) => {
                                    if (value === "filtered" || value === "highlighted") {
                                        setSelectionScope(value);
                                    }
                                }}
                                value={selectionScope}
                            >
                                <FormControlLabel
                                    control={<Radio size="small" />}
                                    label="Filtered rows"
                                    value="filtered"
                                />
                                <FormControlLabel
                                    control={<Radio size="small" />}
                                    label="Highlighted rows"
                                    value="highlighted"
                                />
                            </RadioGroup>
                        </FormControl>
                    </Stack>
                </SectionCard>

                {selectedColumn ? (
                    <TagView
                        columnName={selectedColumn}
                        dataStore={dataStore}
                        selectionScope={selectionScope}
                    />
                ) : (
                    <SectionCard
                        title="Apply Annotations"
                        description="Select or create a multitext column above to start annotating rows."
                    >
                        <Typography color="text.secondary" variant="body2">
                            No annotation column is active yet.
                        </Typography>
                    </SectionCard>
                )}
            </Stack>
        );
    },
);

class AnnotationDialogReact extends BaseDialog {
    root: ReturnType<typeof createMdvPortal>;

    constructor(dataStore: DataStore) {
        super(
            {
                title: `Annotate '${dataStore.name}'`,
                width: 520,
                height: 520,
            },
            null,
        );
        this.root = createMdvPortal(
            <AnnotationDialogComponent dataStore={dataStore} />,
            this.dialog,
            this,
        );
    }

    close(): void {
        super.close();
        this.root.unmount();
    }
}

BaseDialog.experiment["AnnotationDialogReact"] = AnnotationDialogReact;
export default "AnnotationDialogReact loaded";
