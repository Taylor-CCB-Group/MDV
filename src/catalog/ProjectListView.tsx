import {
    Delete,
    DriveFileRenameOutline,
    Image as ImageIcon,
    Info,
    LockPerson,
    MoreVert,
    Upload as UploadIcon,
    Share as ShareIcon
} from "@mui/icons-material";
import {
    IconButton,
    ListItemIcon,
    ListItemText,
    Menu,
    MenuItem,
    Paper,
    Table,
    TableBody,
    TableCell,
    TableContainer,
    TableHead,
    TableRow,
    Typography,
    useTheme,
} from "@mui/material";
import { useCallback, useMemo, useState } from "react";
import ProjectDeleteModal from "./ProjectDeleteModal";
import ProjectInfoModal from "./ProjectInfoModal";
import ProjectRenameModal from "./ProjectRenameModal";
import type { Project } from "./utils/projectUtils";
import type { ProjectAccessType } from "./utils/projectUtils";
import ProjectShareModal from "./ProjectShareModal";
import usePermissions from "./PermissionsContext";

export type Pvoid = Promise<void>;
export type ProjectListViewProps = {
    projects: Project[];
    onDelete: (id: string) => Pvoid;
    onRename: (id: string, name: string) => Pvoid;
    onChangeType: (id: string, type: ProjectAccessType) => Pvoid;
    onExport: (id: string, name: string) => Pvoid;
};
export type MouseEv = React.MouseEvent<HTMLElement>;
const ProjectListView = ({ projects, onDelete, onRename, onExport, onChangeType }: ProjectListViewProps) => {
    const [anchorEl, setAnchorEl] = useState<HTMLElement | null>();
    const [selectedProject, setSelectedProject] = useState<Project>();

    const [isInfoModalOpen, setIsInfoModalOpen] = useState(false);
    const [isRenameModalOpen, setIsRenameModalOpen] = useState(false);
    const [isDeleteModalOpen, setIsDeleteModalOpen] = useState(false);
    const [isShareModalOpen, setIsShareModalOpen] = useState(false);
    const theme = useTheme();
    const { permissions: operationPermissions, isPublicPage } = usePermissions();

    const handleMenuClick = useCallback((event: MouseEv, project: Project) => {
        event.stopPropagation();
        setAnchorEl(event.currentTarget);
        setSelectedProject(project);
    }, []);

    const handleMenuClose = useCallback(() => {
        setAnchorEl(null);
    }, []);

    const handleRowClick = useCallback((projectId: Project['id']) => {
        const base = import.meta.env.DEV ? "http://localhost:5170?dir=/" : "";
        window.location.href = `${base}project/${projectId}`;
    }, []);

    const handleModalOpen = useCallback((modalSetter: React.Dispatch<React.SetStateAction<boolean>>) => {
        handleMenuClose();
        modalSetter(true);
    }, [handleMenuClose]);

    const hasReadme = useCallback((proj: Project) => !!proj.readme, []);

    const hasPermissions = useCallback((proj: Project) => 
        proj.permissions.edit && (operationPermissions.renameProject ||
        operationPermissions.deleteProject ||
        operationPermissions.exportProject ||
        operationPermissions.shareProject)
        , [operationPermissions]);

    // todo: Update this component for more cleaner code
    return (
        <>
            <TableContainer component={Paper}>
                <Table>
                    <TableHead>
                        <TableRow>
                            <TableCell>Name</TableCell>
                            <TableCell>Owner</TableCell>
                            {!isPublicPage && <TableCell>Last Modified</TableCell>}
                            <TableCell>Type</TableCell>
                            <TableCell align="right">Actions</TableCell>
                        </TableRow>
                    </TableHead>
                    <TableBody>
                        {projects.map((project) => (
                            <TableRow
                                key={project.id}
                                hover
                                onClick={() => handleRowClick(project.id)}
                                sx={{ cursor: "pointer" }}
                            >
                                <TableCell>
                                    <div className="flex items-center gap-3">
                                        <ImageIcon
                                            sx={{ color: "text.secondary" }}
                                        />
                                        <Typography>{project.name}</Typography>
                                    </div>
                                </TableCell>
                                <TableCell>{project.owner ? project.owner.join(', ') : ''}</TableCell>
                                {!isPublicPage && <TableCell>{project.lastModified}</TableCell>}
                                <TableCell>{project.type}</TableCell>
                                <TableCell align="right">
                                {(hasReadme(project) || hasPermissions(project)) && 
                                    (
                                        <IconButton
                                            onClick={(e) =>
                                                handleMenuClick(e, project)
                                            }
                                            size="small"
                                        >
                                            <MoreVert />
                                        </IconButton>
                                )}
                                </TableCell>
                            </TableRow>
                        ))}
                    </TableBody>
                </Table>
            </TableContainer>

            {selectedProject && (
                <>
                    <Menu
                        anchorEl={anchorEl}
                        open={Boolean(anchorEl)}
                        onClose={handleMenuClose}
                    >
                        {hasReadme(selectedProject) && (
                            <MenuItem onClick={() => handleModalOpen(setIsInfoModalOpen)}>
                                <ListItemIcon>
                                    <Info fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Project Info</ListItemText>
                            </MenuItem>
                        )}
                        {operationPermissions.renameProject && (
                            <MenuItem 
                                onClick={() => handleModalOpen(setIsRenameModalOpen)}
                                disabled={!selectedProject?.permissions.edit}
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
                                    if (selectedProject) {
                                        handleMenuClose();
                                        onExport(selectedProject.id, selectedProject.name);
                                    }
                                }}
                                disabled={!selectedProject?.permissions.edit}
                            >
                                <ListItemIcon>
                                    <UploadIcon fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Export Project (as *.mdv.zip)</ListItemText>
                            </MenuItem>
                        )}
                        {operationPermissions.shareProject && (
                            <MenuItem 
                                onClick={() => handleModalOpen(setIsShareModalOpen)} 
                                disabled={!selectedProject?.permissions.edit}
                            >
                                <ListItemIcon>
                                    <ShareIcon fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Share Project</ListItemText>
                            </MenuItem>
                        )}
                        {operationPermissions.deleteProject && (
                            <MenuItem 
                                onClick={() => handleModalOpen(setIsDeleteModalOpen)} 
                                sx={{color: theme.palette.error.main}}
                                disabled={!selectedProject?.permissions.edit}
                            >
                                <ListItemIcon>
                                    <Delete color="error" fontSize="small" />
                                </ListItemIcon>
                                <ListItemText>Delete Project</ListItemText>
                            </MenuItem>
                        )}
                    </Menu>

                    {hasReadme(selectedProject) && (
                        <ProjectInfoModal
                        open={isInfoModalOpen}
                        onClose={() => setIsInfoModalOpen(false)}
                        readme={selectedProject.readme}
                    />
                    )}

                    {selectedProject.permissions.edit && (
                        <>
                            {operationPermissions.renameProject && <ProjectRenameModal
                                id={selectedProject.id}
                                name={selectedProject.name}
                                open={isRenameModalOpen}
                                onRename={onRename}
                                onClose={() => setIsRenameModalOpen(false)}
                            />}

                            {operationPermissions.shareProject && <ProjectShareModal
                                open={isShareModalOpen}
                                onClose={() => setIsShareModalOpen(false)}
                                projectId={selectedProject.id}
                            />}

                            {operationPermissions.deleteProject && <ProjectDeleteModal
                                id={selectedProject.id}
                                open={isDeleteModalOpen}
                                onDelete={onDelete}
                                onClose={() => setIsDeleteModalOpen(false)}
                            />}
                        </>
                    )}
                </>
            )}
        </>
    );
};

export default ProjectListView;
