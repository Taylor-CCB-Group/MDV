import {
    Delete as DeleteIcon,
    DriveFileRenameOutline,
    Image as ImageIcon,
    Info,
    LockPerson as LockPersonIcon,
    MoreVert,
    Upload as UploadIcon,
    Settings,
    Share as ShareIcon,
} from "@mui/icons-material";
import {
    Card,
    CardContent,
    CardMedia,
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
import type React from "react";
import { useEffect, useMemo, useState } from "react";
import ProjectDeleteModal from "./ProjectDeleteModal";
import ProjectInfoModal from "./ProjectInfoModal";
import ProjectRenameModal from "./ProjectRenameModal";
import ProjectSettingsModal from "./ProjectSettingsModal";
import type { Permissions, ProjectAccessType } from "./utils/projectUtils";
import ProjectShareModal from "./ProjectShareModal";
import usePermissions from "./PermissionsContext";

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
    onDelete: (id: string) => Promise<void>;
    onRename: (id: string, newName: string) => Promise<void>;
    onChangeType: (
        id: string,
        newType: ProjectAccessType,
    ) => Promise<void>;
    onAddCollaborator: (email: string) => void;
    onExport: (id: string, name: string) => Promise<void>;
    thumbnail?: string;
    readme?: string;
}

const ProjectCard: React.FC<ProjectCardProps> = ({
    id,
    name,
    type,
    lastModified,
    permissions,
    onDelete,
    onRename,
    onChangeType,
    onExport,
    thumbnail,
    readme,
}) => {
    const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
    const [isInfoModalOpen, setIsInfoModalOpen] = useState(false);
    const [isSettingsModalOpen, setIsSettingsModalOpen] = useState(false);
    const [isShareModalOpen, setIsShareModalOpen] = useState(false);
    const [isRenameModalOpen, setIsRenameModalOpen] = useState(false);
    const [isDeleteModalOpen, setIsDeleteModalOpen] = useState(false);
    const [isAccessModalOpen, setIsAccessModalOpen] = useState(false);
    const theme = useTheme();
    const { permissions: operationPermissions, isPublicPage } = usePermissions();
    
    // todo - review how we do stuff like this
    const base = import.meta.env.DEV
        ? "http://localhost:5170?dir=/"
        : "";
    const href = `${base}project/${id}`;
    
    const handleMenuOpen = (event: React.MouseEvent<HTMLButtonElement>) => {
        event.stopPropagation();
        setAnchorEl(event.currentTarget);
    };

    const handleMenuClose = () => {
        setAnchorEl(null);
    };

    const hasReadme = useMemo(() => !!readme, [readme]);
    const hasPermissions = useMemo(
        () =>
            permissions.edit &&
            (operationPermissions.renameProject ||
                operationPermissions.deleteProject ||
                operationPermissions.exportProject ||
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
            }}
            data-testid="project_card"
            data-project-id={id}
        >
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
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "center",
                        bgcolor: thumbnail ? "inherit" : "grey.200",
                    }}
                >
                    {thumbnail ? (
                        <img src={thumbnail} alt="project_thumbnail" style={{ objectFit: "contain", width: "90%", aspectRatio: 3/2 }} />
                    ) : (
                        <ImageIcon sx={{ fontSize: 80, color: "text.secondary" }} />
                    )}
                </CardMedia>

                <Divider />

                <CardContent
                    sx={{
                        padding: "16px",
                        paddingBottom: "12px !important",
                        display: "flex",
                        flexDirection: "column",
                        justifyContent: "flex-end",
                    }}
                >
                    <Tooltip 
                        title={name} 
                        placement="bottom-start"
                        enterDelay={500}
                        slotProps={{
                            popper: {
                                modifiers: [
                                    {
                                        name: 'offset',
                                        options: {
                                            offset: [0, -8],
                                        }
                                    }
                                ]
                            },
                            tooltip: {
                                sx: {
                                    fontSize: "0.9rem",
                                    fontWeight: "normal",
                                }
                            },
                        }}
                    >
                        <Typography
                            gutterBottom
                            variant="h6"
                            component="div"
                            noWrap
                            color="text.primary"
                            sx={{ marginBottom: !isPublicPage ? "4px" : "12px" }}
                        >
                            {name}
                        </Typography>
                    </Tooltip>
                    {!isPublicPage && (
                        <Typography
                            variant="body2"
                            color="text.secondary"
                            sx={{ marginBottom: "0" }}
                        >
                            Last modified: {lastModified}
                        </Typography>
                    )}
                </CardContent>
            </a>

            {(hasReadme || hasPermissions) && 
                (
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
                                bgcolor: "grey.500"
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
                        {operationPermissions.exportProject && (
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
                        )}
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
                                sx={{color: theme.palette.error.main}}
                                disabled={!permissions.edit}
                                data-testid={`project_delete_${id}`}
                            >
                                <ListItemIcon>
                                    <DeleteIcon color="error" fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Delete Project</ListItemText>
                            </MenuItem>
                        )}
                    </Menu>
                </>
            )}

            {hasReadme && (
                <ProjectInfoModal
                    open={isInfoModalOpen}
                    onClose={() => setIsInfoModalOpen(false)}
                    readme={readme}
                />
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

                    {operationPermissions.renameProject && <ProjectRenameModal
                        id={id}
                        name={name}
                        open={isRenameModalOpen}
                        onRename={onRename}
                        onClose={() => setIsRenameModalOpen(false)}
                    />}

                    {operationPermissions.shareProject && <ProjectShareModal
                        open={isShareModalOpen}
                        onClose={() => setIsShareModalOpen(false)}
                        projectId={id}
                    />}

                    {operationPermissions.deleteProject && <ProjectDeleteModal
                        id={id}
                        open={isDeleteModalOpen}
                        onDelete={onDelete}
                        onClose={() => setIsDeleteModalOpen(false)}
                    />}
                </>
            )}
        </Card>
    );
};

export default ProjectCard;
