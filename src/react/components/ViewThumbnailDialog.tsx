import {
    Box,
    Button,
    Card,
    CardActionArea,
    Checkbox,
    Chip,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Divider,
    FormControlLabel,
    IconButton,
    TextField,
    Typography,
    Grid2,
    CircularProgress,
    InputAdornment,
} from "@mui/material";
import {
    Check as CheckIcon,
    Clear as ClearIcon,
    Close as CloseIcon,
    Image as ImageIcon,
} from "@mui/icons-material";
import {
    useCallback,
    useEffect,
    useRef,
    useState,
    type Dispatch,
    type SetStateAction,
} from "react";
import { ErrorBoundary } from "react-error-boundary";
import { stringContainsAll } from "@/lib/utils";
import DebugErrorComponent from "@/charts/dialogs/DebugErrorComponent";

export type ViewThumbnailDialogProps = {
    open: boolean;
    setOpen: Dispatch<SetStateAction<boolean>>;
};

export type ViewImageComponentProps = {
    imgSrc: string;
    viewName: string;
};

export const ViewImageComponent = ({ imgSrc, viewName }: ViewImageComponentProps) => {
    const [hasError, setHasError] = useState(!imgSrc);

    const ZOOM_LEVEL = 1.75;
    const VERTICAL_OFFSET_PERCENT = 25;
    const baseWidth = `${100 / ZOOM_LEVEL}%`;

    return (
        <Box
            sx={{
                width: "100%",
                height: "100%",
                aspectRatio: 3 / 2,
                display: "flex",
                justifyContent: "center",
                alignItems: "flex-start",
                overflow: "hidden",
            }}
        >
            {!hasError ? (
                <img
                    src={imgSrc}
                    alt={`${viewName} snapshot`}
                    style={{
                        width: baseWidth,
                        height: "auto",
                        transform: `scale(${ZOOM_LEVEL}) translateY(${VERTICAL_OFFSET_PERCENT}%)`,
                        transition: "transform 0.2s ease-in-out",
                        transformOrigin: "center top",
                    }}
                    onError={() => setHasError(true)}
                />
            ) : (
                <ImageIcon
                    sx={{
                        fontSize: "5rem",
                        color: "text.secondary",
                    }}
                />
            )}
        </Box>
    );
};

type ViewEntry = { name: string; image: string };

type RenameState = { index: number; value: string } | null;

const ViewThumbnailDialog = ({ open, setOpen }: ViewThumbnailDialogProps) => {
    const viewManager = window.mdv.chartManager.viewManager;
    const config = window.mdv.chartManager.config;
    const isEditMode = config.permission === "edit" || config.permission === "owner";

    const [viewList, setViewList] = useState<ViewEntry[]>([]);
    const [filteredViewList, setFilteredViewList] = useState<ViewEntry[]>([]);
    const [loading, setLoading] = useState(false);
    const [filterText, setFilterText] = useState<string>("");

    // Rename state: which card is being renamed and the current input value
    const [renaming, setRenaming] = useState<RenameState>(null);
    const [renameError, setRenameError] = useState<string>("");
    const [renameSaving, setRenameSaving] = useState(false);

    // Drag-and-drop state
    const dragIndexRef = useRef<number | null>(null);
    const [dragOverIndex, setDragOverIndex] = useState<number | null>(null);
    const [orderChanged, setOrderChanged] = useState(false);
    const [savingOrder, setSavingOrder] = useState(false);

    // Gallery default checkbox
    const [galleryDefault, setGalleryDefault] = useState<boolean>(
        Boolean(config.show_gallery_on_open),
    );

    const getViewList = useCallback(async () => {
        setLoading(true);
        const list = await viewManager.getViewDetails();
        setViewList(list);
        setFilteredViewList(list);
        setLoading(false);
    }, [viewManager]);

    useEffect(() => {
        if (viewManager && open) {
            getViewList();
            setOrderChanged(false);
            setRenaming(null);
            setRenameError("");
        }
    }, [viewManager, getViewList, open]);

    const onClose = () => {
        setOpen(false);
    };

    const handleCardClick = async (name: string) => {
        if (renaming !== null) return;
        viewManager.checkAndChangeView(name);
        onClose();
    };

    const onInputChange = (inputText: string) => {
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

    const startRename = (e: React.MouseEvent, index: number) => {
        e.stopPropagation();
        setRenaming({ index, value: filteredViewList[index].name });
        setRenameError("");
    };

    const cancelRename = () => {
        setRenaming(null);
        setRenameError("");
    };

    const commitRename = async () => {
        if (!renaming) return;
        const trimmed = renaming.value.trim();
        if (!trimmed) {
            setRenameError("Name cannot be empty");
            return;
        }
        const oldName = filteredViewList[renaming.index].name;
        if (trimmed === oldName) {
            cancelRename();
            return;
        }
        if (viewList.some((v) => v.name === trimmed)) {
            setRenameError("A view with that name already exists");
            return;
        }
        setRenameSaving(true);
        try {
            await viewManager.renameView(oldName, trimmed);
            // Update local lists
            const updateEntry = (v: ViewEntry) =>
                v.name === oldName ? { ...v, name: trimmed } : v;
            setViewList((prev) => prev.map(updateEntry));
            setFilteredViewList((prev) => prev.map(updateEntry));
            setRenaming(null);
            setRenameError("");
        } catch {
            setRenameError("Failed to rename view. Please try again.");
        } finally {
            setRenameSaving(false);
        }
    };

    // ── Drag-and-drop helpers ───────────────────────────────────────────────────

    const handleDragStart = (e: React.DragEvent, index: number) => {
        dragIndexRef.current = index;
        e.dataTransfer.effectAllowed = "move";
    };

    const handleDragOver = (e: React.DragEvent, index: number) => {
        e.preventDefault();
        e.dataTransfer.dropEffect = "move";
        setDragOverIndex(index);
    };

    const handleDragLeave = () => {
        setDragOverIndex(null);
    };

    const handleDrop = (e: React.DragEvent, targetIndex: number) => {
        e.preventDefault();
        setDragOverIndex(null);
        const srcIndex = dragIndexRef.current;
        if (srcIndex === null || srcIndex === targetIndex) return;

        const reordered = [...filteredViewList];
        const [moved] = reordered.splice(srcIndex, 1);
        reordered.splice(targetIndex, 0, moved);
        setFilteredViewList(reordered);

        // Also update the master viewList to match
        setViewList(reordered);
        setOrderChanged(true);
        dragIndexRef.current = null;
    };

    const handleDragEnd = () => {
        dragIndexRef.current = null;
        setDragOverIndex(null);
    };

    const handleSaveOrder = async () => {
        setSavingOrder(true);
        try {
            await viewManager.reorderViews(filteredViewList.map((v) => v.name));
            setOrderChanged(false);
        } catch {
            // order stays changed so user can retry
        } finally {
            setSavingOrder(false);
        }
    };

    // ── Gallery default checkbox ────────────────────────────────────────────────

    const handleGalleryDefaultChange = async (checked: boolean) => {
        setGalleryDefault(checked);
        try {
            await viewManager.setGalleryDefault(checked);
        } catch {
            setGalleryDefault(!checked);
        }
    };

    return (
        <Dialog open={open} onClose={onClose} fullWidth maxWidth="xl">
            <DialogTitle>
                Browse View Gallery
                <IconButton
                    aria-label="close"
                    onClick={onClose}
                    sx={{
                        position: "absolute",
                        right: 8,
                        top: 8,
                        color: (theme: any) => theme.palette.grey[500],
                    }}
                >
                    <CloseIcon />
                </IconButton>
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
                                    variant="caption"
                                    color="text.secondary"
                                    sx={{ mb: 2, display: "block" }}
                                >
                                    Drag cards to reorder views. Double-click a name to rename.
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
                                <Grid2 container spacing={4}>
                                    {filteredViewList.map((view, index) => {
                                        const isRenaming =
                                            renaming !== null && renaming.index === index;
                                        const isDragTarget = dragOverIndex === index;

                                        return (
                                            <Grid2
                                                key={`${view.name}-${index}`}
                                                size={{ md: 4, sm: 6 }}
                                                draggable={isEditMode}
                                                onDragStart={(e: React.DragEvent<HTMLDivElement>) => {
                                                    if (isEditMode) handleDragStart(e, index);
                                                }}
                                                onDragOver={(e: React.DragEvent<HTMLDivElement>) => {
                                                    if (isEditMode) handleDragOver(e, index);
                                                }}
                                                onDragLeave={() => {
                                                    if (isEditMode) handleDragLeave();
                                                }}
                                                onDrop={(e: React.DragEvent<HTMLDivElement>) => {
                                                    if (isEditMode) handleDrop(e, index);
                                                }}
                                                onDragEnd={() => {
                                                    if (isEditMode) handleDragEnd();
                                                }}
                                                sx={{
                                                    opacity: isDragTarget ? 0.5 : 1,
                                                    cursor: isEditMode ? "grab" : "default",
                                                    transition: "opacity 0.15s",
                                                }}
                                            >
                                                <Card
                                                    sx={{
                                                        boxShadow: isDragTarget ? 8 : 20,
                                                        bgcolor: "inherit",
                                                        outline: isDragTarget
                                                            ? "2px dashed"
                                                            : "none",
                                                        outlineColor: isDragTarget
                                                            ? "primary.main"
                                                            : "transparent",
                                                    }}
                                                >
                                                    <CardActionArea
                                                        onClick={() =>
                                                            !isRenaming &&
                                                            handleCardClick(view.name)
                                                        }
                                                        sx={{
                                                            cursor: isRenaming
                                                                ? "default"
                                                                : "pointer",
                                                        }}
                                                    >
                                                        <Box
                                                            sx={{
                                                                padding: "1%",
                                                                bgcolor: "background.default",
                                                                position: "relative",
                                                            }}
                                                        >
                                                            {isEditMode && (
                                                                <Chip
                                                                    label={index + 1}
                                                                    size="small"
                                                                    sx={{
                                                                        position: "absolute",
                                                                        top: 8,
                                                                        left: 8,
                                                                        zIndex: 1,
                                                                        fontWeight: "bold",
                                                                        minWidth: 28,
                                                                    }}
                                                                />
                                                            )}
                                                            <ViewImageComponent
                                                                imgSrc={view.image}
                                                                viewName={view.name}
                                                            />
                                                        </Box>
                                                        <Divider />
                                                        <Box
                                                            sx={{
                                                                display: "flex",
                                                                justifyContent: "center",
                                                                alignItems: "center",
                                                                paddingY: "3%",
                                                                paddingX: "4%",
                                                                minHeight: 52,
                                                            }}
                                                            onClick={(e: React.MouseEvent) => {
                                                                if (isEditMode) e.stopPropagation();
                                                            }}
                                                        >
                                                            {isRenaming ? (
                                                                <Box
                                                                    sx={{
                                                                        display: "flex",
                                                                        alignItems: "center",
                                                                        gap: 0.5,
                                                                        width: "100%",
                                                                    }}
                                                                    onClick={(e: React.MouseEvent) =>
                                                                        e.stopPropagation()
                                                                    }
                                                                >
                                                                    <TextField
                                                                        autoFocus
                                                                        size="small"
                                                                        value={renaming.value}
                                                                        error={Boolean(renameError)}
                                                                        helperText={renameError}
                                                                        disabled={renameSaving}
                                                                        onChange={(e: React.ChangeEvent<HTMLInputElement>) =>
                                                                            setRenaming({
                                                                                ...renaming,
                                                                                value: e.target.value,
                                                                            })
                                                                        }
                                                                        onKeyDown={(e: React.KeyboardEvent) => {
                                                                            if (e.key === "Enter")
                                                                                commitRename();
                                                                            if (e.key === "Escape")
                                                                                cancelRename();
                                                                        }}
                                                                        sx={{ flex: 1 }}
                                                                    />
                                                                    <IconButton
                                                                        size="small"
                                                                        color="primary"
                                                                        disabled={renameSaving}
                                                                        onClick={(e: React.MouseEvent) => {
                                                                            e.stopPropagation();
                                                                            commitRename();
                                                                        }}
                                                                        title="Save name"
                                                                    >
                                                                        {renameSaving ? (
                                                                            <CircularProgress
                                                                                size={16}
                                                                            />
                                                                        ) : (
                                                                            <CheckIcon fontSize="small" />
                                                                        )}
                                                                    </IconButton>
                                                                    <IconButton
                                                                        size="small"
                                                                        disabled={renameSaving}
                                                                        onClick={(e: React.MouseEvent) => {
                                                                            e.stopPropagation();
                                                                            cancelRename();
                                                                        }}
                                                                        title="Cancel"
                                                                    >
                                                                        <ClearIcon fontSize="small" />
                                                                    </IconButton>
                                                                </Box>
                                                            ) : (
                                                                <Typography
                                                                    sx={{ fontWeight: "bold" }}
                                                                    onDoubleClick={(e: React.MouseEvent) => {
                                                                        if (isEditMode) startRename(e, index);
                                                                    }}
                                                                    title={
                                                                        isEditMode
                                                                            ? "Double-click to rename"
                                                                            : undefined
                                                                    }
                                                                >
                                                                    {view.name}
                                                                </Typography>
                                                            )}
                                                        </Box>
                                                    </CardActionArea>
                                                </Card>
                                            </Grid2>
                                        );
                                    })}
                                </Grid2>
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
                {isEditMode && orderChanged && (
                    <Button
                        onClick={handleSaveOrder}
                        color="primary"
                        variant="contained"
                        disabled={savingOrder}
                        sx={{ mr: "auto", ml: 2 }}
                    >
                        {savingOrder ? <CircularProgress size={18} sx={{ mr: 1 }} /> : null}
                        Save Order
                    </Button>
                )}
                <Button onClick={onClose} color="primary" sx={{ mr: 2 }}>
                    Cancel
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default ViewThumbnailDialog;
