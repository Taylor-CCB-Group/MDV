import type React from "react";
import { useState } from "react";
import {
    MoreVert,
    Info,
    Settings,
    Share,
    Image as ImageIcon,
} from "@mui/icons-material";
import {
    Card,
    CardContent,
    CardMedia,
    Typography,
    IconButton,
    Menu,
    MenuItem,
    ListItemIcon,
    ListItemText,
} from "@mui/material";
import ProjectInfoModal from "./ProjectInfoModal";
import ProjectSettingsModal from "./ProjectSettingsModal";
import ProjectShareModal from "./ProjectShareModal";

interface ProjectCardProps {
    id: string;
    name: string;
    type: "Editable" | "Read-Only";
    lastModified: string;
    createdAt: string;
    owner: string;
    collaborators: string[];
    numberOfStructures: string;
    numberOfImages: string;
    onDelete: (id: string) => Promise<void>;
    onRename?: (id: string, newName: string) => void;
    onChangeType?: (id: string, newType: "Editable" | "Read-Only") => void;
    onAddCollaborator?: (email: string) => void;
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

    const handleMenuOpen = (event: React.MouseEvent<HTMLButtonElement>) => {
        event.stopPropagation();
        setAnchorEl(event.currentTarget);
    };

    const handleMenuClose = () => {
        setAnchorEl(null);
    };

    const handleCardClick = (event: React.MouseEvent<HTMLDivElement>) => {
        // Prevent navigation if clicking on the menu button or its children
        if (!(event.target as HTMLElement).closest('.menu-button')) {
            window.location.href = `http://localhost:5170?dir=/project/${id}`;
        }
    };

    return (
        <Card
            sx={{
                aspectRatio: "1 / 1",
                width: "100%",
                display: "flex",
                flexDirection: "column",
                position: "relative",
                cursor: "pointer",
            }}
            onClick={handleCardClick}
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
            >
                <MenuItem
                    onClick={(e) => {
                        e.stopPropagation();
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
                    onClick={(e) => {
                        e.stopPropagation();
                        setIsSettingsModalOpen(true);
                        handleMenuClose();
                    }}
                >
                    <ListItemIcon>
                        <Settings fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Project Settings</ListItemText>
                </MenuItem>
                <MenuItem
                    onClick={(e) => {
                        e.stopPropagation();
                        setIsShareModalOpen(true);
                        handleMenuClose();
                    }}
                >
                    <ListItemIcon>
                        <Share fontSize="small" />
                    </ListItemIcon>
                    <ListItemText>Share Project</ListItemText>
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
        </Card>
    );
};

export default ProjectCard;