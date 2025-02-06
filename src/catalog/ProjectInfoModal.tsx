import CloseIcon from "@mui/icons-material/Close";
import {
    Box,
    Dialog,
    DialogContent,
    DialogTitle,
    Divider,
    IconButton,
    List,
    ListItem,
    ListItemText,
    Typography,
} from "@mui/material";
import type React from "react";

export interface ProjectInfoModalProps {
    open: boolean;
    name: string;
    createdAt: string;
    lastModified: string;
    owner: string;
    collaborators: string[];
    numberOfStructures: string;
    numberOfImages: string;
    onClose: () => void;
}

const ProjectInfoModal: React.FC<ProjectInfoModalProps> = ({
    open,
    name,
    createdAt,
    lastModified,
    owner,
    collaborators,
    numberOfImages,
    numberOfStructures,
    onClose,
}) => {
    return (
        <Dialog
            open={open}
            onClose={onClose}
            aria-labelledby="project-info-dialog-title"
            maxWidth="sm"
            fullWidth
        >
            <DialogTitle id="project-info-dialog-title">
                Project Information
                <IconButton
                    aria-label="close"
                    onClick={onClose}
                    sx={{
                        position: "absolute",
                        right: 8,
                        top: 8,
                        color: (theme) => theme.palette.grey[500],
                    }}
                >
                    <CloseIcon />
                </IconButton>
            </DialogTitle>
            <DialogContent dividers>
                <List>
                    <ListItem>
                        <ListItemText primary="Project Name" secondary={name} />
                    </ListItem>
                    <Divider component="li" />
                    <ListItem>
                        <ListItemText
                            primary="Created At"
                            secondary={createdAt}
                        />
                    </ListItem>
                    <Divider component="li" />
                    <ListItem>
                        <ListItemText
                            primary="Last Modified"
                            secondary={lastModified}
                        />
                    </ListItem>
                    <Divider component="li" />
                    <ListItem>
                        <ListItemText
                            primary="Project Owner"
                            secondary={owner}
                        />
                    </ListItem>
                    <Divider component="li" />
                    <ListItem>
                        <ListItemText
                            primary="Number of Structures"
                            secondary={numberOfStructures}
                        />
                    </ListItem>
                    <Divider component="li" />
                    <ListItem>
                        <ListItemText
                            primary="Number of Images"
                            secondary={numberOfImages}
                        />
                    </ListItem>
                    <Divider component="li" />
                    <ListItem>
                        <ListItemText
                            primary="Collaborators"
                            secondary={
                                <Box component="span" sx={{ display: "block" }}>
                                    {collaborators.map((collaborator) => (
                                        <Typography
                                            key={collaborator}
                                            variant="body2"
                                            component="p"
                                        >
                                            {collaborator}
                                        </Typography>
                                    ))}
                                </Box>
                            }
                        />
                    </ListItem>
                </List>
            </DialogContent>
        </Dialog>
    );
};

export default ProjectInfoModal;
