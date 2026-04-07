import {
    useEffect,
    useMemo,
    useState,
    type MouseEvent,
    type PropsWithChildren,
} from "react";
import { observer } from "mobx-react-lite";
import {
    Accordion,
    AccordionContent,
    AccordionItem,
    AccordionTrigger,
} from "@/components/ui/accordion";
import TabHeader from "@/react/components/TabHeader";
import { getNextHighlightedRows } from "@/react/selectionHooks";
import type DataStore from "../../datastore/DataStore.js";
import TagModel, { type TagSelectionScope } from "../../table/TagModel";
import {
    Autocomplete,
    Box,
    Button,
    Checkbox,
    Chip,
    Divider,
    Paper,
    Stack,
    TextField,
    Typography,
} from "@mui/material";
import CheckBoxOutlineBlankIcon from "@mui/icons-material/CheckBoxOutlineBlank";
import CheckBoxIcon from "@mui/icons-material/CheckBox";
import WarningAmberIcon from "@mui/icons-material/WarningAmber";

const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;
const SCOPE_TABS = ["Filtered", "Highlighted"] as const;

type ScopeTab = (typeof SCOPE_TABS)[number];
type AnnotationColumnChoiceState = "empty" | "ok" | "clash" | "readonly";
function isEditableAnnotationColumn(
    column: { datatype?: string; editable?: boolean } | undefined,
) {
    return column?.datatype === "multitext" && column.editable !== false;
}

export function getEditableAnnotationColumnNames(
    columns: Array<{ datatype?: string; editable?: boolean; field?: string; name?: string }>,
) {
    return columns
        .filter(isEditableAnnotationColumn)
        .map((column) => column.field || column.name || "")
        .filter(Boolean)
        .sort((left, right) => left.localeCompare(right));
}

export function getAnnotationColumnChoiceState(
    dataStore: DataStore,
    rawValue: string,
): AnnotationColumnChoiceState {
    const nextValue = rawValue.trim();
    if (!nextValue) {
        return "empty";
    }

    const existingColumn = dataStore.columnIndex[nextValue];
    if (!existingColumn) {
        return "ok";
    }

    if (existingColumn.datatype !== "multitext") {
        return "clash";
    }

    return existingColumn.editable === false ? "readonly" : "ok";
}

function getDefaultAnnotationColumn(dataStore: DataStore) {
    const editableColumnNames = getEditableAnnotationColumnNames(dataStore.columns);
    if (editableColumnNames.includes("__tags")) {
        return "__tags";
    }
    return editableColumnNames[0] || "";
}

function getSelectionScope(tab: ScopeTab): TagSelectionScope {
    return tab === "Highlighted" ? "highlighted" : "filtered";
}

function getChoiceStateHelperText(choiceState: AnnotationColumnChoiceState) {
    if (choiceState === "clash") {
        return "A non-multitext column already uses that name.";
    }
    if (choiceState === "readonly") {
        return "That multitext column is read-only and can't be edited here.";
    }
    if (choiceState === "empty") {
        return "Enter a column name to select or create an editable annotation column.";
    }
    return "Names become new editable multitext annotation columns when needed.";
}

function useTagModelState(
    dataStore: DataStore,
    columnName: string,
    selectionScope: TagSelectionScope,
) {
    const [tagModel, setTagModel] = useState<TagModel | null>(null);
    const [loadError, setLoadError] = useState<string | null>(null);
    const [, setRevision] = useState(0);

    useEffect(() => {
        let cancelled = false;
        let activeModel: TagModel | null = null;
        let listener: (() => void) | null = null;

        setTagModel(null);
        setLoadError(null);

        if (!columnName) {
            return;
        }

        TagModel.create(dataStore, columnName, selectionScope)
            .then((model) => {
                if (cancelled) {
                    model.dispose();
                    return;
                }

                activeModel = model;
                listener = () => {
                    setRevision((n) => n + 1);
                };

                model.addListener(listener);
                setTagModel(model);
                setRevision((n) => n + 1);
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
            if (activeModel && listener) {
                activeModel.removeListener(listener);
            }
            activeModel?.dispose();
        };
    }, [columnName, dataStore, selectionScope]);

    return {
        tagModel,
        loadError,
    };
}

function SetupSection({
    dataStore,
    columnInput,
    onColumnInputChange,
    onApplyColumnChoice,
    selectedColumn,
}: {
    dataStore: DataStore;
    columnInput: string;
    onColumnInputChange: (value: string) => void;
    onApplyColumnChoice: (value: string) => void;
    selectedColumn: string;
}) {
    const [isOpen, setIsOpen] = useState(() => !selectedColumn);
    const editableColumnNames = getEditableAnnotationColumnNames(dataStore.columns);
    const editableColumnNamesKey = editableColumnNames.join("\0");
    const choiceState = getAnnotationColumnChoiceState(dataStore, columnInput);

    useEffect(() => {
        if (!selectedColumn) {
            setIsOpen(true);
        }
    }, [selectedColumn]);

    return (
        <Paper
            className="w-full"
            variant="outlined"
            sx={{
                width: "100%",
                backgroundColor: (theme) =>
                    theme.palette.mode === "dark"
                        ? "rgba(255,255,255,0.03)"
                        : "rgba(0,0,0,0.02)",
                borderColor: (theme) => theme.palette.divider,
                borderWidth: 1,
                borderStyle: "solid",
                borderRadius: 1,
                overflow: "hidden",
                paddingX: 2,
                paddingBottom: 0,
            }}
        >
            <Accordion
                type="single"
                collapsible
                className="w-full"
                value={isOpen ? "annotation-setup" : ""}
                onValueChange={(value) => {
                    setIsOpen(value === "annotation-setup");
                }}
            >
                <AccordionItem value="annotation-setup" className="border-b-0">
                    <AccordionTrigger className="w-full py-3 hover:no-underline">
                        <Box sx={{ textAlign: "left" }}>
                            <Typography variant="h6">Annotation Setup</Typography>
                            <Typography color="text.secondary" variant="body2">
                                Choose or create an editable multitext column for tag
                                annotations.
                            </Typography>
                        </Box>
                    </AccordionTrigger>
                    <AccordionContent className="pb-5 pt-2">
                        <Stack spacing={2}>
                            <Stack direction={{ sm: "row", xs: "column" }} spacing={1}>
                                <Autocomplete
                                    freeSolo
                                    fullWidth
                                    size="small"
                                    options={editableColumnNames}
                                    inputValue={columnInput}
                                    onChange={(_, value) => {
                                        if (typeof value === "string") {
                                            onColumnInputChange(value);
                                            if (
                                                getAnnotationColumnChoiceState(dataStore, value) ===
                                                "ok"
                                            ) {
                                                onApplyColumnChoice(value);
                                                setIsOpen(false);
                                            }
                                        }
                                    }}
                                    onInputChange={(_, value) => {
                                        onColumnInputChange(value);
                                    }}
                                    renderInput={(params) => (
                                        <TextField
                                            {...params}
                                            error={choiceState === "clash" || choiceState === "readonly"}
                                            helperText={getChoiceStateHelperText(choiceState)}
                                            label="Column name"
                                            onKeyDown={(event) => {
                                                if (event.key !== "Enter") {
                                                    return;
                                                }
                                                event.preventDefault();
                                                if (choiceState === "ok") {
                                                    onApplyColumnChoice(columnInput);
                                                    setIsOpen(false);
                                                }
                                            }}
                                        />
                                    )}
                                />
                                <Button
                                    disabled={choiceState !== "ok"}
                                    onClick={() => {
                                        onApplyColumnChoice(columnInput);
                                        setIsOpen(false);
                                    }}
                                    sx={{ minWidth: 112 }}
                                    variant="contained"
                                >
                                    Use Column
                                </Button>
                            </Stack>

                            <Typography color="text.secondary" variant="body2">
                                Editable annotation columns:{" "}
                                {editableColumnNamesKey
                                    ? editableColumnNames.join(", ")
                                    : "none yet"}
                            </Typography>
                        </Stack>
                    </AccordionContent>
                </AccordionItem>
            </Accordion>
        </Paper>
    );
}

function WorkspaceCard({
    title,
    description,
    children,
}: PropsWithChildren<{ title: string; description?: string }>) {
    return (
        <Paper
            className="w-full"
            variant="outlined"
            sx={{
                width: "100%",
                borderColor: "divider",
                overflow: "hidden",
            }}
        >
            <Box sx={{ px: 2.5, pt: 2.5, pb: 2 }}>
                <Typography variant="h6">{title}</Typography>
                {description ? (
                    <Typography color="text.secondary" variant="body2">
                        {description}
                    </Typography>
                ) : null}
            </Box>
            <Divider />
            {children}
        </Paper>
    );
}

function TagEditorPanel({
    dataStore,
    columnName,
    selectionScope,
    tagModel,
    loadError,
}: {
    dataStore: DataStore;
    columnName: string;
    selectionScope: TagSelectionScope;
    tagModel: TagModel | null;
    loadError: string | null;
}) {
    const [draftTag, setDraftTag] = useState("");
    const [tagError, setTagError] = useState<string | null>(null);

    const { tagList, tagsInSelection, selectionCount } = tagModel
        ? tagModel.getAnnotationViewState()
        : {
              tagList: new Set<string>(),
              tagsInSelection: new Set<string>(),
              selectionCount: 0,
          };
    const { showFilteredScopeCoversWholeTableWarning, showNoSelectionWarning } =
        tagModel
            ? tagModel.getAnnotationWarningFlags()
            : {
                  showFilteredScopeCoversWholeTableWarning: false,
                  showNoSelectionWarning: true,
              };

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

    useEffect(() => {
        setDraftTag("");
        setTagError(null);
    }, [columnName, selectionScope]);

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

    const selectionLabel =
        selectionScope === "highlighted" ? "highlighted" : "filtered";
    const showFilteredWarning = showFilteredScopeCoversWholeTableWarning;
    const showNoRowsWarning = showNoSelectionWarning;
    const noRowsLabel =
        selectionScope === "highlighted"
            ? "No highlighted rows match right now."
            : "No filtered rows match right now.";
    const noRowsHelperText =
        selectionScope === "highlighted"
            ? "Highlight some rows before applying tags in this scope."
            : "Adjust the current filters, or switch to Highlighted to work from a focused subset.";

    const addTag = (rawTag: string) => {
        const nextTag = rawTag.trim();
        if (!nextTag) {
            return;
        }
        try {
            tagModel.setTag(nextTag, true);
            setTagError(null);
            setDraftTag("");
        } catch (error) {
            setTagError(error instanceof Error ? error.message : String(error));
        }
    };

    const updateTag = (tag: string, tagValue: boolean) => {
        try {
            tagModel.setTag(tag, tagValue);
            setTagError(null);
        } catch (error) {
            setTagError(error instanceof Error ? error.message : String(error));
        }
    };

    const highlightMatchingRows = (
        tag: string,
        event: MouseEvent<HTMLButtonElement>,
    ) => {
        event.stopPropagation();
        const matchingRows = tagModel.getMatchingRowIndices(tag);
        const currentHighlights = dataStore.getHighlightedData()?.slice() || [];
        dataStore.dataHighlighted(
            getNextHighlightedRows(matchingRows, currentHighlights, event),
            tagModel,
        );
    };

    return (
        <Stack spacing={2.5}>
            <Box>
                <Typography variant="subtitle1">
                    Annotating {selectionCount} {selectionLabel} row
                    {selectionCount === 1 ? "" : "s"}
                </Typography>
                <Typography color="text.secondary" variant="body2">
                    Column: {columnName}
                </Typography>
            </Box>

            {showNoRowsWarning ? (
                <Paper
                    sx={{
                        alignItems: "flex-start",
                        backgroundColor: (theme) =>
                            theme.palette.mode === "dark"
                                ? "rgba(245, 158, 11, 0.12)"
                                : "rgba(245, 158, 11, 0.08)",
                        borderColor: (theme) =>
                            theme.palette.mode === "dark"
                                ? "rgba(245, 158, 11, 0.35)"
                                : "rgba(217, 119, 6, 0.28)",
                        display: "flex",
                        gap: 1.25,
                        p: 1.5,
                    }}
                    variant="outlined"
                >
                    <WarningAmberIcon color="warning" fontSize="small" sx={{ mt: 0.2 }} />
                    <Box>
                        <Typography variant="body2">{noRowsLabel}</Typography>
                        <Typography color="text.secondary" variant="caption">
                            {noRowsHelperText}
                        </Typography>
                    </Box>
                </Paper>
            ) : null}

            {showFilteredWarning ? (
                <Paper
                    sx={{
                        alignItems: "flex-start",
                        backgroundColor: (theme) =>
                            theme.palette.mode === "dark"
                                ? "rgba(245, 158, 11, 0.12)"
                                : "rgba(245, 158, 11, 0.08)",
                        borderColor: (theme) =>
                            theme.palette.mode === "dark"
                                ? "rgba(245, 158, 11, 0.35)"
                                : "rgba(217, 119, 6, 0.28)",
                        display: "flex",
                        gap: 1.25,
                        p: 1.5,
                    }}
                    variant="outlined"
                >
                    <WarningAmberIcon color="warning" fontSize="small" sx={{ mt: 0.2 }} />
                    <Box>
                        <Typography variant="body2">
                            Filtered mode will apply to every row right now.
                        </Typography>
                        <Typography color="text.secondary" variant="caption">
                            Use Highlighted for targeted edits when you only mean to tag a subset.
                        </Typography>
                    </Box>
                </Paper>
            ) : null}

            <Stack direction={{ sm: "row", xs: "column" }} spacing={1}>
                <Autocomplete
                    freeSolo
                    fullWidth
                    size="small"
                    options={availableTags}
                    inputValue={draftTag}
                    onInputChange={(_, value) => {
                        setDraftTag(value);
                        if (tagError) {
                            setTagError(null);
                        }
                    }}
                    onChange={(_, value) => {
                        if (typeof value === "string") {
                            addTag(value);
                        }
                    }}
                    renderInput={(params) => (
                        <TextField
                            {...params}
                            error={Boolean(tagError)}
                            helperText={tagError || " "}
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
                    )}
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
                                    onDelete={() => updateTag(tag, false)}
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
                <Typography color="text.secondary" sx={{ mb: 1.5 }} variant="body2">
                    Click a row to add or remove a tag on the current scope. Use
                    Highlight to replace the current highlight set, hold Shift to add,
                    or Ctrl/Cmd to toggle matching rows.
                </Typography>
                {availableTags.length === 0 ? (
                    <Typography color="text.secondary" variant="body2">
                        Add the first annotation item for this column above.
                    </Typography>
                ) : (
                    <Paper
                        sx={{
                            maxHeight: 260,
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
                                        updateTag(
                                            tag,
                                            indeterminate ? true : !entireSelection,
                                        )
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
                                        sx={{ p: 0.5, pointerEvents: "none" }}
                                    />
                                    <Box sx={{ flexGrow: 1, minWidth: 0 }}>
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
                                    <Button
                                        onClick={(event) => highlightMatchingRows(tag, event)}
                                        size="small"
                                        variant="text"
                                    >
                                        Highlight
                                    </Button>
                                </Box>
                            );
                        })}
                    </Paper>
                )}
            </Box>
        </Stack>
    );
}

function AnnotationWorkspace({
    dataStore,
    selectedColumn,
}: {
    dataStore: DataStore;
    selectedColumn: string;
}) {
    const [activeTab, setActiveTab] = useState<ScopeTab>(() =>
        (dataStore.getHighlightedData()?.length || 0) > 0 ? "Highlighted" : "Filtered",
    );
    const filteredState = useTagModelState(dataStore, selectedColumn, "filtered");
    const highlightedState = useTagModelState(
        dataStore,
        selectedColumn,
        "highlighted",
    );

    if (!selectedColumn) {
        return (
            <WorkspaceCard
                title="Apply Annotations"
                description="Choose or create an editable multitext column above to start annotating rows."
            >
                <Box sx={{ px: 2.5, py: 2.5 }}>
                    <Typography color="text.secondary" variant="body2">
                        No annotation column is active yet.
                    </Typography>
                </Box>
            </WorkspaceCard>
        );
    }

    const selectionScope = getSelectionScope(activeTab);
    const activeState =
        selectionScope === "highlighted" ? highlightedState : filteredState;

    return (
        <WorkspaceCard title="Apply Annotations">
            <Box sx={{ px: 2.5, pt: 1 }}>
                <TabHeader
                    activeTab={activeTab}
                    setActiveTab={setActiveTab}
                    tabs={[...SCOPE_TABS]}
                />
            </Box>
            <Box sx={{ px: 2.5, py: 2.5 }}>
                <TagEditorPanel
                    columnName={selectedColumn}
                    dataStore={dataStore}
                    loadError={activeState.loadError}
                    selectionScope={selectionScope}
                    tagModel={activeState.tagModel}
                />
            </Box>
        </WorkspaceCard>
    );
}

export const AnnotationDialogComponent = observer(
    ({ dataStore }: { dataStore: DataStore }) => {
        const [selectedColumn, setSelectedColumn] = useState(() =>
            getDefaultAnnotationColumn(dataStore),
        );
        const [columnInput, setColumnInput] = useState(selectedColumn);

        const editableColumnNames = getEditableAnnotationColumnNames(dataStore.columns);
        const editableColumnNamesKey = editableColumnNames.join("\0");

        useEffect(() => {
            if (
                selectedColumn &&
                getAnnotationColumnChoiceState(dataStore, selectedColumn) === "ok"
            ) {
                return;
            }

            const fallbackColumn = getDefaultAnnotationColumn(dataStore);
            if (fallbackColumn !== selectedColumn) {
                setSelectedColumn(fallbackColumn);
            }
            if (!columnInput.trim() || columnInput === selectedColumn) {
                setColumnInput(fallbackColumn);
            }
        }, [columnInput, dataStore, editableColumnNamesKey, selectedColumn]);

        const applyColumnChoice = (rawValue: string) => {
            const nextValue = rawValue.trim();
            if (getAnnotationColumnChoiceState(dataStore, nextValue) !== "ok") {
                return;
            }
            setSelectedColumn(nextValue);
            setColumnInput(nextValue);
        };

        return (
            <Stack spacing={2}>
                <SetupSection
                    columnInput={columnInput}
                    dataStore={dataStore}
                    onApplyColumnChoice={applyColumnChoice}
                    onColumnInputChange={setColumnInput}
                    selectedColumn={selectedColumn}
                />
                <AnnotationWorkspace
                    dataStore={dataStore}
                    selectedColumn={selectedColumn}
                />
            </Stack>
        );
    },
);
