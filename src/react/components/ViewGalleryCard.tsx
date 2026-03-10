import {
    Box,
    Card,
    CardActionArea,
    Chip,
    Divider,
    Grid2,
    IconButton,
    TextField,
    Typography,
    CircularProgress,
} from "@mui/material";
import { Check as CheckIcon, Clear as ClearIcon } from "@mui/icons-material";
import { useSortable } from "@dnd-kit/sortable";
import { CSS } from "@dnd-kit/utilities";
import ViewImageComponent from "./ViewImageComponent";

export type ViewGalleryCardView = { name: string; image: string };

export type ViewGalleryCardProps = {
    view: ViewGalleryCardView;
    index: number;
    isEditMode: boolean;
    canReorder: boolean;
    isRenaming: boolean;
    renameValue: string;
    renameError: string;
    renameSaving: boolean;
    onCardClick: (name: string) => void;
    onStartRename: (e: React.MouseEvent, index: number) => void;
    onRenameChange: (value: string) => void;
    onRenameKeyDown: (e: React.KeyboardEvent) => void;
    onCommitRename: () => void;
    onCancelRename: () => void;
};

export function ViewGalleryCard({
    view,
    index,
    isEditMode,
    canReorder,
    isRenaming,
    renameValue,
    renameError,
    renameSaving,
    onCardClick,
    onStartRename,
    onRenameChange,
    onRenameKeyDown,
    onCommitRename,
    onCancelRename,
}: ViewGalleryCardProps) {
    const {
        attributes,
        listeners,
        setNodeRef,
        transform,
        transition,
        isDragging,
    } = useSortable({
        id: view.name,
        // Disabled when filter is on or not in edit mode.
        disabled: !canReorder,
    });

    // Apply transform and transition so the card moves smoothly during drag.
    const style = {
        transform: CSS.Translate.toString(transform),
        transition,
    };

    return (
        <Grid2
            ref={setNodeRef}
            size={{ md: 4, sm: 6 }}
            sx={{ opacity: isDragging ? 0.5 : 1 }}
            {...(canReorder ? { ...attributes, ...listeners } : {})}
            style={style}
        >
            <Card
                sx={{
                    boxShadow: 20,
                    bgcolor: "inherit",
                }}
            >
                <CardActionArea
                    onClick={() => !isRenaming && onCardClick(view.name)}
                    sx={{ cursor: isRenaming ? "default" : "pointer", display: "block" }}
                >
                    <Box sx={{ padding: "1%", bgcolor: "background.default", position: "relative" }}>
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
                                    fontSize: "0.95rem",
                                }}
                            />
                        )}
                        <ViewImageComponent imgSrc={view.image} viewName={view.name} />
                    </Box>
                </CardActionArea>
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
                >
                    {isRenaming ? (
                        <Box
                            sx={{ display: "flex", alignItems: "center", gap: 0.5, width: "100%" }}
                            onClick={(e) => e.stopPropagation()}
                        >
                            <TextField
                                autoFocus
                                size="small"
                                value={renameValue}
                                error={Boolean(renameError)}
                                helperText={renameError}
                                disabled={renameSaving}
                                onChange={(e) => onRenameChange(e.target.value)}
                                onKeyDown={onRenameKeyDown}
                                sx={{ flex: 1 }}
                            />
                            <IconButton
                                size="small"
                                color="primary"
                                disabled={renameSaving}
                                onClick={(e) => {
                                    e.stopPropagation();
                                    onCommitRename();
                                }}
                                aria-label="Save name"
                                title="Save name"
                            >
                                {renameSaving ? (
                                    <CircularProgress size={16} />
                                ) : (
                                    <CheckIcon fontSize="small" color="success" />
                                )}
                            </IconButton>
                            <IconButton
                                size="small"
                                disabled={renameSaving}
                                onClick={(e) => {
                                    e.stopPropagation();
                                    onCancelRename();
                                }}
                                aria-label="Cancel rename"
                                title="Cancel"
                            >
                                <ClearIcon fontSize="small" color="error" />
                            </IconButton>
                        </Box>
                    ) : (
                        <Typography
                            component={"span"}
                            sx={{ fontWeight: "bold" }}
                            onDoubleClick={(e) => isEditMode && onStartRename(e, index)}
                            title={isEditMode ? "Double-click to rename" : undefined}
                            tabIndex={isEditMode ? 0 : undefined} // make label focusable
                            role={isEditMode ? "button" : undefined}
                            onKeyDown={(e) => {
                                if (!isEditMode) return;
                                if (e.key === "Enter") {
                                    e.preventDefault();
                                    onStartRename(e as unknown as React.MouseEvent, index);
                                }
                            }}
                        >
                            {view.name}
                        </Typography>
                    )}
                </Box>
            </Card>
        </Grid2>
    );
}
