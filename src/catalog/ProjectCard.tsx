import {
    Delete as DeleteIcon,
    DriveFileRenameOutline,
    Image as ImageIcon,
    Info,
    LockPerson as LockPersonIcon,
    MoreVert,
    Upload as UploadIcon,
    Settings,
    Share,
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
    Typography,
    useTheme,
} from "@mui/material";
import type React from "react";
import { useEffect, useState } from "react";
import ProjectDeleteModal from "./ProjectDeleteModal";
import ProjectInfoModal from "./ProjectInfoModal";
import ProjectRenameModal from "./ProjectRenameModal";
import ProjectSettingsModal from "./ProjectSettingsModal";
import ProjectShareModal from "./ProjectShareModal";
import type { ProjectAccessType } from "./utils/projectUtils";
import axios from "axios";

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
    onExport: (id: string, name: string) => Promise<void>;
    thumbnail?: string;
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
    onExport,
    thumbnail,
}) => {
    const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
    const [isInfoModalOpen, setIsInfoModalOpen] = useState(false);
    const [isSettingsModalOpen, setIsSettingsModalOpen] = useState(false);
    const [isShareModalOpen, setIsShareModalOpen] = useState(false);
    const [isRenameModalOpen, setIsRenameModalOpen] = useState(false);
    const [isDeleteModalOpen, setIsDeleteModalOpen] = useState(false);
    const [isAccessModalOpen, setIsAccessModalOpen] = useState(false);
    const theme = useTheme();
    
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
            <a
                href={href}
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
                    <Typography
                        gutterBottom
                        variant="h6"
                        component="div"
                        noWrap
                        color="text.primary"
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
            </a>

            <IconButton
                className="menu-button"
                aria-label="project options"
                onClick={handleMenuOpen}
                sx={{
                    position: "absolute",
                    top: 8,
                    right: 8,
                    bgcolor: "inherit",
                    "&:hover": { bgcolor: "inherit" },
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
                        handleMenuClose();
                        onExport(id, name)
                    }}
                >
                    <ListItemIcon>
                        <UploadIcon fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Export Project (as *.zip)</ListItemText>
                </MenuItem>
                <MenuItem
                    onClick={() => {
                        setIsDeleteModalOpen(true);
                        handleMenuClose();
                    }}
                    sx={{color: theme.palette.error.main}}
                >
                    <ListItemIcon>
                        <DeleteIcon color="error" fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Delete Project</ListItemText>
                </MenuItem>
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
        </Card>
    );
};

export default ProjectCard;
