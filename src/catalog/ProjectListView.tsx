import {
    Delete,
    DriveFileRenameOutline,
    Image as ImageIcon,
    Info,
    LockPerson,
    MoreVert,
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
} from "@mui/material";
import { useState } from "react";
import ProjectAccessModal from "./ProjectAccessModal";
import ProjectDeleteModal from "./ProjectDeleteModal";
import ProjectInfoModal from "./ProjectInfoModal";
import ProjectRenameModal from "./ProjectRenameModal";
import type { Project } from "./utils/projectUtils";
import type { ProjectAccessType } from "./utils/projectUtils";

type Pvoid = Promise<void>;
type ProjectListViewProps = {
    projects: Project[];
    onDelete: (id: string) => Pvoid;
    onRename: (id: string, name: string) => Pvoid;
    onChangeType: (id: string, type: ProjectAccessType) => Pvoid;
};
const ProjectListView = ({ projects, onDelete, onRename, onChangeType }: ProjectListViewProps) => {
    const [anchorEl, setAnchorEl] = useState<HTMLAnchorElement | null>();
    const [selectedProject, setSelectedProject] = useState<Project>();

    const [isInfoModalOpen, setIsInfoModalOpen] = useState(false);
    const [isRenameModalOpen, setIsRenameModalOpen] = useState(false);
    const [isDeleteModalOpen, setIsDeleteModalOpen] = useState(false);
    const [isAccessModalOpen, setIsAccessModalOpen] = useState(false);

    const handleMenuClick = (event, project) => {
        event.stopPropagation();
        setAnchorEl(event.currentTarget);
        setSelectedProject(project);
    };

    const handleMenuClose = () => {
        setAnchorEl(null);
    };

    const handleRowClick = (projectId) => {
        const base = import.meta.env.DEV ? "http://localhost:5170?dir=/" : "";
        window.location.href = `${base}project/${projectId}`;
    };

    const handleModalOpen = (modalSetter) => {
        handleMenuClose();
        modalSetter(true);
    };

    return (
        <>
            <TableContainer component={Paper}>
                <Table>
                    <TableHead>
                        <TableRow>
                            <TableCell>Name</TableCell>
                            <TableCell>Owner</TableCell>
                            <TableCell>Last Modified</TableCell>
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
                                <TableCell>{project.owner}</TableCell>
                                <TableCell>{project.lastModified}</TableCell>
                                <TableCell>{project.type}</TableCell>
                                <TableCell align="right">
                                    <IconButton
                                        onClick={(e) =>
                                            handleMenuClick(e, project)
                                        }
                                        size="small"
                                    >
                                        <MoreVert />
                                    </IconButton>
                                </TableCell>
                            </TableRow>
                        ))}
                    </TableBody>
                </Table>
            </TableContainer>

            <Menu
                anchorEl={anchorEl}
                open={Boolean(anchorEl)}
                onClose={handleMenuClose}
            >
                <MenuItem onClick={() => handleModalOpen(setIsInfoModalOpen)}>
                    <ListItemIcon>
                        <Info fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Project Info</ListItemText>
                </MenuItem>
                <MenuItem onClick={() => handleModalOpen(setIsRenameModalOpen)}>
                    <ListItemIcon>
                        <DriveFileRenameOutline fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Rename Project</ListItemText>
                </MenuItem>
                <MenuItem onClick={() => handleModalOpen(setIsDeleteModalOpen)}>
                    <ListItemIcon>
                        <Delete fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Delete Project</ListItemText>
                </MenuItem>
                <MenuItem onClick={() => handleModalOpen(setIsAccessModalOpen)}>
                    <ListItemIcon>
                        <LockPerson fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Project Access</ListItemText>
                </MenuItem>
            </Menu>

            {selectedProject && (
                <>
                    <ProjectInfoModal
                        open={isInfoModalOpen}
                        onClose={() => setIsInfoModalOpen(false)}
                        name={selectedProject.name}
                        createdAt={selectedProject.createdAt}
                        lastModified={selectedProject.lastModified}
                        owner={selectedProject.owner}
                        collaborators={selectedProject.collaborators}
                        numberOfStructures={selectedProject.numberOfStructures}
                        numberOfImages={selectedProject.numberOfImages}
                    />

                    <ProjectRenameModal
                        id={selectedProject.id}
                        name={selectedProject.name}
                        open={isRenameModalOpen}
                        onRename={onRename}
                        onClose={() => setIsRenameModalOpen(false)}
                    />

                    <ProjectDeleteModal
                        id={selectedProject.id}
                        open={isDeleteModalOpen}
                        onDelete={onDelete}
                        onClose={() => setIsDeleteModalOpen(false)}
                    />

                    <ProjectAccessModal
                        id={selectedProject.id}
                        type={selectedProject.type}
                        open={isAccessModalOpen}
                        onChangeType={onChangeType}
                        onClose={() => setIsAccessModalOpen(false)}
                    />
                </>
            )}
        </>
    );
};

export default ProjectListView;
