import {
    Close as CloseIcon,
    ContentCopy as ContentCopyIcon,
    Email as EmailIcon,
} from "@mui/icons-material";
import {
    Box,
    Button,
    Dialog,
    DialogContent,
    DialogTitle,
    IconButton,
    Snackbar,
    TextField,
    Typography,
} from "@mui/material";
import type React from "react";
import { useState } from "react";

interface ProjectShareModalProps {
    open: boolean;
    onClose: () => void;
    onAddCollaborator: (email: string) => void;
    projectId: string;
}

const ProjectShareModal: React.FC<ProjectShareModalProps> = ({
    open,
    onClose,
    onAddCollaborator,
    projectId,
}) => {
    const [email, setEmail] = useState("");
    const [copySuccess, setCopySuccess] = useState(false);

    const projectUrl = `${window.location.origin}/projects/${projectId}`;

    const handleAdd = () => {
        if (email) {
            onAddCollaborator(email);
            setEmail(""); // Clear input field after adding
        }
    };

    const handleCopyLink = () => {
        navigator.clipboard.writeText(projectUrl).then(() => {
            setCopySuccess(true);
        });
    };

    return (
        <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
            <DialogTitle>
                Share Project
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
            <DialogContent>
                <Box sx={{ mt: 2 }}>
                    <Typography variant="subtitle1" gutterBottom>
                        Invite Collaborator
                    </Typography>
                    <Box sx={{ display: "flex", gap: 1, mb: 3 }}>
                        <TextField
                            fullWidth
                            size="small"
                            type="email"
                            value={email}
                            onChange={(e) => setEmail(e.target.value)}
                            placeholder="example@mail.com"
                            variant="outlined"
                        />
                        <Button
                            variant="contained"
                            color="primary"
                            onClick={handleAdd}
                            startIcon={<EmailIcon />}
                        >
                            Invite
                        </Button>
                    </Box>
                    <Typography variant="subtitle1" gutterBottom>
                        Deep link
                    </Typography>
                    <Box sx={{ display: "flex", gap: 1 }}>
                        <TextField
                            fullWidth
                            size="small"
                            value={projectUrl}
                            InputProps={{
                                readOnly: true,
                            }}
                            variant="outlined"
                            onClick={(e) =>
                                (e.target as HTMLInputElement).select()
                            }
                        />
                        <Button
                            variant="outlined"
                            onClick={handleCopyLink}
                            startIcon={<ContentCopyIcon />}
                        >
                            Copy
                        </Button>
                    </Box>
                </Box>
            </DialogContent>
            <Snackbar
                open={copySuccess}
                autoHideDuration={2000}
                onClose={() => setCopySuccess(false)}
                message="Link copied to clipboard"
            />
        </Dialog>
    );
};

export default ProjectShareModal;
