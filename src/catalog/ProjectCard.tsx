import type React from "react";
import { useState } from "react";
import {
    MoreVert,
    Info,
    Settings,
    Share,
    Edit,
    Delete,
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
    Box,
} from "@mui/material";
import ProjectInfoModal from "./ProjectInfoModal";
import ProjectSettingsModal from "./ProjectSettingsModal";
import ProjectShareModal from "./ProjectShareModal";

interface ProjectCardProps {
    id: string;
    name: string;
    type: "Editable" | "Read-Only";
    lastModified?: string;
    imageUrl?: string;
    createdAt: string;
    owner: string;
    collaborators: string[];
    numberOfStructures: string;
    numberOfImages: string;
    onDelete?: (id: string) => void;
    onRename?: (id: string, newName: string) => void;
    onChangeType?: (id: string, newType: "Editable" | "Read-Only") => void;
    onAddCollaborator?: (email: string) => void;
}

const ProjectCard: React.FC<ProjectCardProps> = ({
    id,
    name,
    type,
    lastModified = "Not available",
    imageUrl,
    createdAt,
    owner,
    collaborators,
    numberOfStructures,
    numberOfImages,
    onDelete = () => {},
    onRename = () => {},
    onChangeType = () => {},
    onAddCollaborator = () => {},
}) => {
    const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);
    const [isInfoModalOpen, setIsInfoModalOpen] = useState(false);
    const [isSettingsModalOpen, setIsSettingsModalOpen] = useState(false);
    const [isShareModalOpen, setIsShareModalOpen] = useState(false);

    const handleMenuOpen = (event: React.MouseEvent<HTMLButtonElement>) => {
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
            {/* Options Button */}
            <IconButton
                aria-label="project options"
                onClick={handleMenuOpen}
                sx={{
                    position: "absolute",
                    top: 8,
                    right: 8,
                    bgcolor: "background.paper",
                    "&:hover": { bgcolor: "action.hover" },
                }}
            >
                <MoreVert />
            </IconButton>

            {/* Image */}
            <CardMedia
                component="div"
                sx={{
                    aspectRatio: "1 / 1", // Ensures the media section is also square
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                    bgcolor: "grey.200",
                }}
            >
                {imageUrl ? (
                    <img
                        src={imageUrl}
                        alt={name}
                        style={{
                            maxHeight: "100%",
                            maxWidth: "100%",
                            objectFit: "contain",
                        }}
                    />
                ) : (
                    <ImageIcon sx={{ fontSize: 80, color: "text.secondary" }} />
                )}
            </CardMedia>

            {/* Project Details */}
            <CardContent
                sx={{
                    padding: "16px", // Default padding
                    paddingBottom: "12px", // Override the bottom padding applied by MUI
                    display: "flex",
                    flexDirection: "column",
                    justifyContent: "flex-end",
                    "&:last-child": {
                        paddingBottom: "12px", // Ensure no padding at the bottom for the last child
                    },
                }}
            >
                <Typography
                    gutterBottom
                    variant="h6"
                    component="div"
                    noWrap
                    sx={{ marginBottom: "4px" }} // Adjust the margin between the text and the bottom
                >
                    {name}
                </Typography>
                <Typography
                    variant="body2"
                    color="text.secondary"
                    sx={{ marginBottom: "0" }} // Remove bottom margin to reduce space
                >
                    Last modified: {lastModified}
                </Typography>
            </CardContent>

            {/* Menu */}
            <Menu
                anchorEl={anchorEl}
                open={Boolean(anchorEl)}
                onClose={handleMenuClose}
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
                    onClick={() => {
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

            {/* Modals */}
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
