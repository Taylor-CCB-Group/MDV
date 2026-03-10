import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    FormControlLabel,
    IconButton,
    TextField,
    Typography,
    Grid2,
    CircularProgress,
    InputAdornment,
    Checkbox,
} from "@mui/material";
import { Clear as ClearIcon, Close as CloseIcon } from "@mui/icons-material";
import { useCallback, useEffect, useState, type Dispatch, type SetStateAction } from "react";
import { ErrorBoundary } from "react-error-boundary";
import {
    DndContext,
    closestCenter,
    PointerSensor,
    useSensor,
    useSensors,
    type DragEndEvent,
} from "@dnd-kit/core";
import { arrayMove, SortableContext, rectSortingStrategy } from "@dnd-kit/sortable";
import { stringContainsAll } from "@/lib/utils";
import DebugErrorComponent from "@/charts/dialogs/DebugErrorComponent";
import { useChartManager } from "../hooks";
import { ViewGalleryCard } from "./ViewGalleryCard";
import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";

export type ViewThumbnailDialogProps = {
    open: boolean;
    setOpen: Dispatch<SetStateAction<boolean>>;
};

type ViewEntry = { name: string; image: string };

type RenameState = { name: string; value: string } | null;

const ViewThumbnailDialog = ({ open, setOpen }: ViewThumbnailDialogProps) => {
    const cm = useChartManager();
    const viewManager = cm.viewManager;
    const config = cm.config;
    const isEditMode = config.permission === "edit" || config.permission === "owner";

    const [viewList, setViewList] = useState<ViewEntry[]>([]);
    const [filteredViewList, setFilteredViewList] = useState<ViewEntry[]>([]);
    const [loading, setLoading] = useState(false);
    const [filterText, setFilterText] = useState<string>("");

    // Rename state: which card is being renamed and the current input value
    const [renaming, setRenaming] = useState<RenameState>(null);
    const [renameError, setRenameError] = useState<string>("");
    const [renameSaving, setRenameSaving] = useState(false);

    const [orderChanged, setOrderChanged] = useState(false);
    const [savingOrder, setSavingOrder] = useState(false);

    // Gallery default checkbox
    const [galleryDefault, setGalleryDefault] = useState<boolean>(
        Boolean(config.show_gallery_on_open),
    );

    const isFilterActive = filterText.trim() !== "";
    const canReorder = isEditMode && !isFilterActive;

    const getViewList = useCallback(async () => {
        setLoading(true);
        // Fetch the view list
        const list = await viewManager.getViewDetails();
        setViewList(list);
        setFilteredViewList(list);
        setLoading(false);
    }, [viewManager]);

    useEffect(() => {
        // When dialog opens, refresh the list and reset edit/drag state.
        if (viewManager && open) {
            getViewList();
            setOrderChanged(false);
            setRenaming(null);
            setRenameError("");
        }
    }, [viewManager, getViewList, open]);

    const onClose = () => {
        cancelRename();
        setFilterText("");
        setOpen(false);
    };

    // Switch to the selected view. Ignore click if user is in the middle of renaming.
    const handleCardClick = async (name: string) => {
        if (renaming !== null) return;
        viewManager.checkAndChangeView(name);
        onClose();
    };

    // Filter views by searched text
    const onInputChange = (inputText: string) => {
        setRenaming(null);
        setFilterText(inputText);
        const input = inputText.toLowerCase().split(" ");
        const tempList = viewList.filter((view) => {
            const name = view.name.toLowerCase();
            return stringContainsAll(input, name);
        });
        setFilteredViewList(tempList);
    };

    const handleClearInput = () => {
        setFilterText("");
        setFilteredViewList(viewList);
    };

    // ── Rename helpers ──────────────────────────────────────────────────────────

    const startRename = (e: React.MouseEvent, viewName: string) => {
        e.stopPropagation();
        setRenaming({ name: viewName, value: viewName });
        setRenameError("");
    };

    const cancelRename = () => {
        setRenaming(null);
        setRenameError("");
    };

    // Validate the input then update both viewList and filteredViewList.
    const commitRename = async () => {
        if (!renaming) return;
        const trimmed = renaming.value.trim();
        if (!trimmed) {
            setRenameError("Name cannot be empty");
            return;
        }
        const oldName = renaming.name;

        // same name
        if (trimmed === oldName) {
            cancelRename();
            return;
        }

        // name already exists
        if (viewList.some((v) => v.name === trimmed)) {
            setRenameError("A view with that name already exists");
            return;
        }
        setRenameSaving(true);

        try {
            await viewManager.renameView(oldName, trimmed);
            // Update the local states
            const updateEntry = (v: ViewEntry) =>
                v.name === oldName ? { ...v, name: trimmed } : v;
            setViewList((prev) => prev.map(updateEntry));
            setFilteredViewList((prev) => prev.map(updateEntry));
            // Reset the rename states
            setRenaming(null);
            setRenameError("");
        } catch {
            setRenameError("Failed to rename view. Please try again.");
        } finally {
            setRenameSaving(false);
        }
    };

    // Only start a drag after the pointer moves 8px; avoids triggering drag on image click.
    const sensors = useSensors(
        useSensor(PointerSensor, { activationConstraint: { distance: 8 } }),
    );

    // On drop: resolve indices by view name (id), reorder with arrayMove, update both lists and mark order changed.
    const handleDragEnd = useCallback(
        (event: DragEndEvent) => {
            const { active, over } = event;
            if (!over || active.id === over.id || !canReorder) return;
            const oldIndex = filteredViewList.findIndex((v) => v.name === active.id);
            const newIndex = filteredViewList.findIndex((v) => v.name === over.id);
            if (oldIndex === -1 || newIndex === -1) return;
            const reordered = arrayMove(filteredViewList, oldIndex, newIndex);
            setFilteredViewList(reordered);
            setViewList(reordered);
            setOrderChanged(true);
        },
        [filteredViewList, canReorder],
    );

    // Persist the current order to the backend. Uses full viewList so order is always complete.
    const handleSaveOrder = async () => {
        setSavingOrder(true);
        try {
            await viewManager.reorderViews(viewList.map((v) => v.name));
            setOrderChanged(false);
        } catch {
            // order stays changed so user can retry
        } finally {
            setSavingOrder(false);
        }
    };

    const handleResetOrder = async () => {
        await getViewList();
        setOrderChanged(false);
    };

    // ── Gallery default checkbox ────────────────────────────────────────────────

    const handleGalleryDefaultChange = async (checked: boolean) => {
        setGalleryDefault(checked);
        try {
            await viewManager.setGalleryDefault(checked);
        } catch {
            console.error("Failed to change gallery default");
            setGalleryDefault(!checked);
        }
    };

    return (
        <Dialog open={open} onClose={onClose} fullWidth maxWidth="xl">
            <DialogTitle>
                Browse View Gallery
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers sx={{ height: "95vh" }}>
                <ErrorBoundary
                    FallbackComponent={({ error }: { error: Error }) => (
                        <DebugErrorComponent
                            error={error}
                            title="Error displaying view gallery"
                        />
                    )}
                >
                    {loading ? (
                        <Box
                            sx={{
                                display: "flex",
                                justifyContent: "center",
                                alignItems: "center",
                                height: "100%",
                                width: "100%",
                            }}
                        >
                            <CircularProgress />
                        </Box>
                    ) : (
                        <Box sx={{ flexGrow: 1 }}>
                            <Box sx={{ mt: 1, mb: 3 }}>
                                <TextField
                                    fullWidth
                                    placeholder="Enter the view name to filter the views"
                                    value={filterText}
                                    onChange={(e: React.ChangeEvent<HTMLInputElement>) => onInputChange(e.target.value)}
                                    slotProps={{
                                        input: {
                                            endAdornment: (
                                                <InputAdornment position="end">
                                                    <IconButton onClick={handleClearInput}>
                                                        <ClearIcon />
                                                    </IconButton>
                                                </InputAdornment>
                                            ),
                                        },
                                    }}
                                />
                            </Box>

                            {isEditMode && (
                                <Typography
                                    variant="body1"
                                    color="textDisabled"
                                    sx={{ mb: 3, display: "flex", justifyContent: "center" }}
                                >
                                    {canReorder
                                        ? "Drag cards to reorder views. Double-click a name to rename."
                                        : "Clear the filter to reorder views. Double-click a name to rename."}
                                </Typography>
                            )}

                            {filteredViewList.length === 0 ? (
                                <Box
                                    sx={{
                                        display: "flex",
                                        justifyContent: "center",
                                        alignItems: "center",
                                        width: "100%",
                                    }}
                                >
                                    <Typography variant="h6">No views found</Typography>
                                </Box>
                            ) : (
                                <DndContext
                                    sensors={sensors}
                                    collisionDetection={closestCenter}
                                    onDragEnd={handleDragEnd}
                                >
                                    <SortableContext
                                        items={filteredViewList.map((v) => v.name)}
                                        strategy={rectSortingStrategy}
                                        disabled={!canReorder}
                                    >
                                        <Grid2 container spacing={4}>
                                            {filteredViewList.map((view, index) => (
                                                <ViewGalleryCard
                                                    key={view.name}
                                                    view={view}
                                                    index={index}
                                                    isEditMode={isEditMode}
                                                    canReorder={canReorder}
                                                    isRenaming={renaming !== null && renaming.name === view.name}
                                                    renameValue={renaming?.name === view.name ? renaming.value : ""}
                                                    renameError={renameError}
                                                    renameSaving={renameSaving}
                                                    onCardClick={handleCardClick}
                                                    onStartRename={startRename}
                                                    onRenameChange={(value) =>
                                                        setRenaming((r) => (r ? { ...r, value } : null))
                                                    }
                                                    onRenameKeyDown={(e) => {
                                                        if (e.key === "Enter") commitRename();
                                                        if (e.key === "Escape") {
                                                            e.stopPropagation();
                                                            e.preventDefault();
                                                            cancelRename();
                                                        }
                                                    }}
                                                    onCommitRename={commitRename}
                                                    onCancelRename={cancelRename}
                                                />
                                            ))}
                                        </Grid2>
                                    </SortableContext>
                                </DndContext>
                            )}

                            {isEditMode && (
                                <Box sx={{ mt: 3, pt: 2, borderTop: 1, borderColor: "divider" }}>
                                    <FormControlLabel
                                        control={
                                            <Checkbox
                                                checked={galleryDefault}
                                                onChange={(e: React.ChangeEvent<HTMLInputElement>) =>
                                                    handleGalleryDefaultChange(e.target.checked)
                                                }
                                                size="small"
                                            />
                                        }
                                        label="Show Gallery View by default when opening this project"
                                    />
                                </Box>
                            )}
                        </Box>
                    )}
                </ErrorBoundary>
            </DialogContent>
            <DialogActions>
                {orderChanged && canReorder && (
                    <Box sx={{ display: "flex", justifyContent: "flex-start", mr: "auto" }}>
                        <Button
                            onClick={handleSaveOrder}
                            color="primary"
                            variant="contained"
                            disabled={savingOrder}
                            sx={{ ml: 2 }}
                        >
                            {savingOrder ? <CircularProgress size={18} sx={{ mr: 1 }} /> : null}
                            Save Order
                        </Button>
                        <Button color="error" sx={{ ml: 1 }} onClick={handleResetOrder}>
                            Reset
                        </Button>
                    </Box>
                )}
                <Button onClick={onClose} color="error" sx={{ mr: 2 }}>
                    Close
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default ViewThumbnailDialog;
