import {
    Delete as DeleteIcon,
    DriveFileRenameOutline,
    Image as ImageIcon,
    Info,
    MoreVert,
    Share as ShareIcon,
    WarningAmber,
} from "@mui/icons-material";
import {
    Box,
    Card,
    CardContent,
    CardMedia,
    Checkbox,
    Chip,
    Divider,
    IconButton,
    ListItemIcon,
    ListItemText,
    Menu,
    MenuItem,
    Tooltip,
    Typography,
    useTheme,
} from "@mui/material";
import { keyframes } from "@mui/material/styles";
import type React from "react";
import { useMemo, useRef, useState } from "react";
import { buildProjectUrl } from "@/utils/mdvRouting";
import usePermissions from "./PermissionsContext";
import ProjectDeleteModal from "./ProjectDeleteModal";
import ProjectInfoModal from "./ProjectInfoModal";
import ProjectRenameModal from "./ProjectRenameModal";
import ProjectSettingsModal from "./ProjectSettingsModal";
import ProjectShareModal from "./ProjectShareModal";
import type { Permissions, ProjectAccessType } from "./utils/projectUtils";

// The full title grows upward from its resting position. The negative insets leave room
// for the box-shadow, so nothing visibly changes when the clip-path is removed afterwards.
const revealTitle = keyframes`
    from { clip-path: inset(100% 0 0 0); opacity: 0; }
    to { clip-path: inset(-12px -12px 0 -12px); opacity: 1; }
`;
const hideTitle = keyframes`
    from { clip-path: inset(-12px -12px 0 -12px); opacity: 1; }
    to { clip-path: inset(100% 0 0 0); opacity: 0; }
`;

const prefersReducedMotion = () =>
    typeof window !== "undefined" && window.matchMedia?.("(prefers-reduced-motion: reduce)").matches;

const titleSx = {
    fontSize: "1.05rem",
    lineHeight: 1.3,
    overflowWrap: "anywhere",
} as const;

/** Long titles wrap to three lines, then truncate; hovering the card shows the whole name. */
const clampedTitleSx = {
    overflow: "hidden",
    display: "-webkit-box",
    WebkitBoxOrient: "vertical",
    WebkitLineClamp: 3,
} as const;

export interface ProjectCardProps {
    id: string;
    name: string;
    type: ProjectAccessType;
    lastModified: string;
    createdAt: string;
    owner: string[];
    collaborators: string[];
    numberOfStructures: string;
    numberOfImages: string;
    permissions: Permissions;
    writable: boolean;
    onDelete: (id: string) => Promise<void>;
    onRename: (id: string, newName: string) => Promise<void>;
    onChangeType: (id: string, newType: ProjectAccessType) => Promise<void>;
    onAddCollaborator: (email: string) => void;
    onExport: (id: string, name: string) => Promise<void>;
    thumbnail?: string;
    readme?: string;
    selectionMode?: boolean;
    selected?: boolean;
    onSelectionChange?: (id: string, selected: boolean) => void;
}

const ProjectCard: React.FC<ProjectCardProps> = ({
    id,
    name,
    type,
    permissions,
    writable,
    onDelete,
    onRename,
    onChangeType,
    onExport,
    thumbnail,
    readme,
    selectionMode = false,
    selected = false,
    onSelectionChange,
}) => {
    const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
    const [isInfoModalOpen, setIsInfoModalOpen] = useState(false);
    const [isSettingsModalOpen, setIsSettingsModalOpen] = useState(false);
    const [isShareModalOpen, setIsShareModalOpen] = useState(false);
    const [isRenameModalOpen, setIsRenameModalOpen] = useState(false);
    const [isDeleteModalOpen, setIsDeleteModalOpen] = useState(false);
    const [isAccessModalOpen, setIsAccessModalOpen] = useState(false);
    // "closing" keeps the overlay mounted while its exit animation plays.
    const [titlePhase, setTitlePhase] = useState<"closed" | "open" | "closing">("closed");
    const titleRef = useRef<HTMLDivElement>(null);
    const theme = useTheme();
    const { permissions: operationPermissions, isPublicPage } = usePermissions();

    // todo - review how we do stuff like this
    const href = buildProjectUrl(id);

    const handleMenuOpen = (event: React.MouseEvent<HTMLButtonElement>) => {
        event.stopPropagation();
        setAnchorEl(event.currentTarget);
    };

    // Only expand when the clamped title is actually cut off.
    const expandTitleIfTruncated = () => {
        const el = titleRef.current;
        if (el && el.scrollHeight > el.clientHeight + 1) setTitlePhase("open");
    };

    const collapseTitle = () => {
        setTitlePhase((phase) => (phase === "open" && !prefersReducedMotion() ? "closing" : "closed"));
    };

    const handleMenuClose = () => {
        setAnchorEl(null);
    };

    const showReadOnlyChip = !isPublicPage && !writable;

    const hasReadme = useMemo(() => !!readme, [readme]);
    const hasPermissions = useMemo(
        () =>
            permissions.edit &&
            (operationPermissions.renameProject ||
                operationPermissions.deleteProject ||
                // operationPermissions.exportProject ||
                operationPermissions.shareProject),
        [permissions.edit, operationPermissions],
    );

    // todo: Update this component for more cleaner code
    return (
        <Card
            sx={{
                aspectRatio: "1 / 1",
                width: "100%",
                display: "flex",
                flexDirection: "column",
                position: "relative",
                transition: "transform 150ms ease-out, box-shadow 150ms ease-out",
                "&:hover, &:has(a:focus-visible)": {
                    transform: "translateY(-2px)",
                    boxShadow: 6,
                },
                // Keep the shadow change as a cue, but drop the movement and animation.
                "@media (prefers-reduced-motion: reduce)": {
                    transition: "none",
                    "&:hover, &:has(a:focus-visible)": {
                        transform: "none",
                    },
                },
            }}
            onMouseEnter={expandTitleIfTruncated}
            onMouseLeave={collapseTitle}
            onFocus={expandTitleIfTruncated}
            onBlur={(event) => {
                if (!event.currentTarget.contains(event.relatedTarget)) collapseTitle();
            }}
            data-testid="project_card"
            data-project-id={id}
        >
            {selectionMode &&
                operationPermissions.deleteProject &&
                permissions.edit &&
                permissions.owner &&
                onSelectionChange && (
                    <Checkbox
                        checked={selected}
                        onChange={(event) => onSelectionChange(id, event.target.checked)}
                        onClick={(event) => event.stopPropagation()}
                        inputProps={{ "aria-label": `select project ${name}` }}
                        sx={{
                            position: "absolute",
                            top: 8,
                            left: 8,
                            bgcolor: "background.paper",
                            zIndex: 1,
                        }}
                    />
                )}
            <a
                href={href}
                style={{
                    cursor: "pointer",
                    flexGrow: 1,
                    display: "flex",
                    flexDirection: "column",
                }}
                aria-label={`open project ${name} (${id})`}
                data-testid={`project_open_${id}`}
            >
                <CardMedia
                    component="div"
                    sx={{
                        flexGrow: 1,
                        minHeight: 0,
                        overflow: "hidden",
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "center",
                        bgcolor: thumbnail ? "inherit" : "grey.200",
                    }}
                >
                    {thumbnail ? (
                        <img
                            src={thumbnail}
                            alt="project_thumbnail"
                            style={{ objectFit: "contain", width: "90%", maxHeight: "100%", aspectRatio: 3 / 2 }}
                        />
                    ) : (
                        <ImageIcon sx={{ fontSize: 80, color: "text.secondary" }} />
                    )}
                </CardMedia>

                <Divider />

                <CardContent
                    sx={{
                        padding: "16px",
                        paddingBottom: "16px !important",
                        display: "flex",
                        flexDirection: "column",
                        justifyContent: "flex-end",
                    }}
                >
                    <Box sx={{ position: "relative", marginBottom: showReadOnlyChip ? "8px" : 0 }}>
                        <Typography
                            ref={titleRef}
                            variant="h6"
                            component="div"
                            color="text.primary"
                            sx={{ ...titleSx, ...clampedTitleSx }}
                        >
                            {name}
                        </Typography>
                        {titlePhase !== "closed" && (
                            // Full name, drawn over the thumbnail so the card layout doesn't shift.
                            <Typography
                                aria-hidden
                                variant="h6"
                                component="div"
                                color="text.primary"
                                sx={{
                                    ...titleSx,
                                    position: "absolute",
                                    bottom: 0,
                                    left: "-16px",
                                    right: "-16px",
                                    px: "16px",
                                    pt: 1,
                                    maxHeight: 220,
                                    overflowY: "auto",
                                    zIndex: 1,
                                    bgcolor: "background.paper",
                                    backgroundImage: "var(--Paper-overlay)",
                                    boxShadow: "0 -4px 8px -4px rgba(0, 0, 0, 0.3)",
                                    animation:
                                        titlePhase === "closing"
                                            ? `${hideTitle} 120ms ease-in forwards`
                                            : `${revealTitle} 160ms ease-out`,
                                    "@media (prefers-reduced-motion: reduce)": { animation: "none" },
                                }}
                                onAnimationEnd={() =>
                                    setTitlePhase((phase) => (phase === "closing" ? "closed" : phase))
                                }
                            >
                                {name}
                            </Typography>
                        )}
                    </Box>
                    {showReadOnlyChip && (
                        <Tooltip title="The MDV server cannot write this project's files. It will open read-only.">
                            <Chip
                                icon={<WarningAmber />}
                                label="Server read-only"
                                color="warning"
                                size="small"
                                sx={{ alignSelf: "flex-start" }}
                            />
                        </Tooltip>
                    )}
                </CardContent>
            </a>

            {(hasReadme || hasPermissions) && (
                <>
                    <IconButton
                        className="menu-button"
                        aria-label="project options"
                        onClick={handleMenuOpen}
                        sx={{
                            position: "absolute",
                            top: 8,
                            right: 8,
                            bgcolor: "background.paper",
                            "&:hover": {
                                bgcolor: "grey.500",
                            },
                            zIndex: 1,
                        }}
                        data-testid={`project_menu_${id}`}
                    >
                        <MoreVert />
                    </IconButton>

                    <Menu
                        anchorEl={anchorEl}
                        open={Boolean(anchorEl)}
                        onClose={handleMenuClose}
                        onClick={(e) => e.stopPropagation()}
                    >
                        {hasReadme && (
                            <MenuItem
                                onClick={() => {
                                    setIsInfoModalOpen(true);
                                    handleMenuClose();
                                }}
                                data-testid={`project_info_${id}`}
                            >
                                <ListItemIcon>
                                    <Info fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Project Info</ListItemText>
                            </MenuItem>
                        )}
                        {operationPermissions.renameProject && (
                            <MenuItem
                                onClick={() => {
                                    setIsRenameModalOpen(true);
                                    handleMenuClose();
                                }}
                                disabled={!permissions.edit}
                                data-testid={`project_rename_${id}`}
                            >
                                <ListItemIcon>
                                    <DriveFileRenameOutline fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Rename Project</ListItemText>
                            </MenuItem>
                        )}
                        {/* {operationPermissions.exportProject && (
                            <MenuItem
                                onClick={() => {
                                    handleMenuClose();
                                    onExport(id, name)
                                }}
                                disabled={!permissions.edit}
                                data-testid={`project_export_${id}`}
                            >
                                <ListItemIcon>
                                    <UploadIcon fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Export Project (as *.mdv.zip)</ListItemText>
                            </MenuItem>
                        )} */}
                        {operationPermissions.shareProject && (
                            <MenuItem
                                onClick={() => {
                                    setIsShareModalOpen(true);
                                    handleMenuClose();
                                }}
                                disabled={!permissions.edit}
                                data-testid={`project_share_${id}`}
                            >
                                <ListItemIcon>
                                    <ShareIcon fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Share Project</ListItemText>
                            </MenuItem>
                        )}
                        {operationPermissions.deleteProject && (
                            <MenuItem
                                onClick={() => {
                                    setIsDeleteModalOpen(true);
                                    handleMenuClose();
                                }}
                                sx={{ color: theme.palette.error.main }}
                                disabled={!permissions.edit}
                                data-testid={`project_delete_${id}`}
                            >
                                <ListItemIcon>
                                    <DeleteIcon color="error" fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Move to Recycle Bin</ListItemText>
                            </MenuItem>
                        )}
                    </Menu>
                </>
            )}

            {hasReadme && (
                <ProjectInfoModal open={isInfoModalOpen} onClose={() => setIsInfoModalOpen(false)} readme={readme} />
            )}

            {permissions.edit && (
                <>
                    <ProjectSettingsModal
                        id={id}
                        name={name}
                        type={type}
                        open={isSettingsModalOpen}
                        onRename={onRename}
                        onChangeType={onChangeType}
                        onDelete={onDelete}
                        onClose={() => setIsSettingsModalOpen(false)}
                    />

                    {operationPermissions.renameProject && (
                        <ProjectRenameModal
                            id={id}
                            name={name}
                            open={isRenameModalOpen}
                            onRename={onRename}
                            onClose={() => setIsRenameModalOpen(false)}
                        />
                    )}

                    {operationPermissions.shareProject && isShareModalOpen && (
                        <ProjectShareModal
                            open={isShareModalOpen}
                            onClose={() => setIsShareModalOpen(false)}
                            projectId={id}
                        />
                    )}

                    {operationPermissions.deleteProject && (
                        <ProjectDeleteModal
                            id={id}
                            open={isDeleteModalOpen}
                            onDelete={onDelete}
                            onClose={() => setIsDeleteModalOpen(false)}
                        />
                    )}
                </>
            )}
        </Card>
    );
};

export default ProjectCard;
