import {
    Delete as DeleteIcon,
    DriveFileRenameOutline,
    Image as ImageIcon,
    Info,
    LockPerson as LockPersonIcon,
    MoreVert,
    Settings,
    Share,
} from "@mui/icons-material";
import {
    Card,
    CardContent,
    CardMedia,
    IconButton,
    ListItemIcon,
    ListItemText,
    Menu,
    MenuItem,
    Typography,
} from "@mui/material";
import type React from "react";
import { useEffect, useState } from "react";
import ProjectAccessModal from "./ProjectAccessModal";
import ProjectDeleteModal from "./ProjectDeleteModal";
import ProjectInfoModal from "./ProjectInfoModal";
import ProjectRenameModal from "./ProjectRenameModal";
import ProjectSettingsModal from "./ProjectSettingsModal";
import ProjectShareModal from "./ProjectShareModal";
import type { ProjectAccessType } from "./utils/projectUtils";

export interface ProjectCardProps {
    id: string;
    name: string;
    type: ProjectAccessType;
    lastModified: string;
    createdAt: string;
    owner: string;
    collaborators: string[];
    numberOfStructures: string;
    numberOfImages: string;
    onDelete: (id: string) => Promise<void>;
    onRename: (id: string, newName: string) => Promise<void>;
    onChangeType: (
        id: string,
        newType: ProjectAccessType,
    ) => Promise<void>;
    onAddCollaborator: (email: string) => void;
}

const ProjectCard: React.FC<ProjectCardProps> = ({
    id,
    name,
    type,
    lastModified,
    createdAt,
    owner,
    collaborators,
    numberOfStructures,
    numberOfImages,
    onDelete,
    onRename,
    onChangeType,
    onAddCollaborator,
}) => {
    const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
    const [isInfoModalOpen, setIsInfoModalOpen] = useState(false);
    const [isSettingsModalOpen, setIsSettingsModalOpen] = useState(false);
    const [isShareModalOpen, setIsShareModalOpen] = useState(false);
    const [isRenameModalOpen, setIsRenameModalOpen] = useState(false);
    const [isDeleteModalOpen, setIsDeleteModalOpen] = useState(false);
    const [isAccessModalOpen, setIsAccessModalOpen] = useState(false);
    const [shouldNavigate, setShouldNavigate] = useState(false);

    useEffect(() => {
        if (shouldNavigate) {
            // todo - review how we do stuff like this
            const base = import.meta.env.DEV
                ? "http://localhost:5170?dir=/"
                : "";
            window.location.href = `${base}project/${id}`;
        }
    }, [shouldNavigate, id]);

    const handleMenuOpen = (event: React.MouseEvent<HTMLButtonElement>) => {
        event.stopPropagation();
        setAnchorEl(event.currentTarget);
    };

    const handleMenuClose = () => {
        setAnchorEl(null);
    };

    const handleCardClick = (event: React.MouseEvent<HTMLDivElement>) => {
        setShouldNavigate(true);
    };

    return (
        <Card
            sx={{
                aspectRatio: "1 / 1",
                width: "100%",
                display: "flex",
                flexDirection: "column",
                position: "relative",
            }}
        >
            <div
                onClick={handleCardClick}
                style={{
                    cursor: "pointer",
                    flexGrow: 1,
                    display: "flex",
                    flexDirection: "column",
                }}
            >
                <CardMedia
                    component="div"
                    sx={{
                        flexGrow: 1,
                        display: "flex",
                        alignItems: "center",
                        justifyContent: "center",
                        bgcolor: "grey.200",
                    }}
                >
                    <ImageIcon sx={{ fontSize: 80, color: "text.secondary" }} />
                </CardMedia>

                <CardContent
                    sx={{
                        padding: "16px",
                        paddingBottom: "12px !important",
                        display: "flex",
                        flexDirection: "column",
                        justifyContent: "flex-end",
                    }}
                >
                    <Typography
                        gutterBottom
                        variant="h6"
                        component="div"
                        noWrap
                        sx={{ marginBottom: "4px" }}
                    >
                        {name}
                    </Typography>
                    <Typography
                        variant="body2"
                        color="text.secondary"
                        sx={{ marginBottom: "0" }}
                    >
                        Last modified: {lastModified}
                    </Typography>
                </CardContent>
            </div>

            <IconButton
                className="menu-button"
                aria-label="project options"
                onClick={handleMenuOpen}
                sx={{
                    position: "absolute",
                    top: 8,
                    right: 8,
                    bgcolor: "background.paper",
                    "&:hover": { bgcolor: "action.hover" },
                    zIndex: 1,
                }}
            >
                <MoreVert />
            </IconButton>

            <Menu
                anchorEl={anchorEl}
                open={Boolean(anchorEl)}
                onClose={handleMenuClose}
                onClick={(e) => e.stopPropagation()}
            >
                <MenuItem
                    onClick={() => {
                        setIsInfoModalOpen(true);
                        handleMenuClose();
                    }}
                >
                    <ListItemIcon>
                        <Info fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Project Info</ListItemText>
                </MenuItem>
                <MenuItem
                    onClick={() => {
                        setIsRenameModalOpen(true);
                        handleMenuClose();
                    }}
                >
                    <ListItemIcon>
                        <DriveFileRenameOutline fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Rename Project</ListItemText>
                </MenuItem>
                <MenuItem
                    onClick={() => {
                        setIsDeleteModalOpen(true);
                        handleMenuClose();
                    }}
                >
                    <ListItemIcon>
                        <DeleteIcon fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Delete Project</ListItemText>
                </MenuItem>
                <MenuItem
                    onClick={() => {
                        setIsAccessModalOpen(true);
                        handleMenuClose();
                    }}
                >
                    <ListItemIcon>
                        <LockPersonIcon fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Project Access</ListItemText>
                </MenuItem>
                {/* <MenuItem
                    onClick={() => {
                        setIsSettingsModalOpen(true);
                        handleMenuClose();
                    }}
                >
                    <ListItemIcon>
                        <Settings fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Project Settings</ListItemText>
                </MenuItem> */}
                {/* <MenuItem
                    onClick={() => {
                        setIsShareModalOpen(true);
                        handleMenuClose();
                    }}
                >
                    <ListItemIcon>
                        <Share fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Share Project</ListItemText>
                </MenuItem> */}
            </Menu>

            <ProjectInfoModal
                open={isInfoModalOpen}
                onClose={() => setIsInfoModalOpen(false)}
                name={name}
                createdAt={createdAt}
                lastModified={lastModified}
                owner={owner}
                collaborators={collaborators}
                numberOfStructures={numberOfStructures}
                numberOfImages={numberOfImages}
            />

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

            <ProjectShareModal
                open={isShareModalOpen}
                onClose={() => setIsShareModalOpen(false)}
                onAddCollaborator={onAddCollaborator}
                projectId={id}
            />

            <ProjectRenameModal
                id={id}
                name={name}
                open={isRenameModalOpen}
                onRename={onRename}
                onClose={() => setIsRenameModalOpen(false)}
            />

            <ProjectDeleteModal
                id={id}
                open={isDeleteModalOpen}
                onDelete={onDelete}
                onClose={() => setIsDeleteModalOpen(false)}
            />

            <ProjectAccessModal
                id={id}
                type={type}
                open={isAccessModalOpen}
                onChangeType={onChangeType}
                onClose={() => setIsAccessModalOpen(false)}
            />
        </Card>
    );
};

export default ProjectCard;
