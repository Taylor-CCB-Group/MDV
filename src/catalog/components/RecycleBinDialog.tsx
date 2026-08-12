import DeleteForeverIcon from "@mui/icons-material/DeleteForever";
import DeleteOutlineIcon from "@mui/icons-material/DeleteOutline";
import FolderOffOutlinedIcon from "@mui/icons-material/FolderOffOutlined";
import {
    alpha,
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Divider,
    Paper,
    Stack,
    Typography,
    useTheme,
} from "@mui/material";
import { useEffect, useState } from "react";
import { DialogCloseIconButton } from "../ProjectRenameModal";
import type { RecycledProject } from "../utils/projectUtils";

interface RecycleBinDialogProps {
    open: boolean;
    onClose: () => void;
    fetchProjects: () => Promise<RecycledProject[]>;
    onDeleteAll: () => Promise<{
        deletedProjectIds: number[];
        failures: { projectId: number; message: string }[];
    }>;
    onRestore: (projectId: string) => Promise<{ restoredProjectId: number }>;
}

const formatDeletedAt = (deletedTimestamp?: string) => {
    if (!deletedTimestamp) {
        return "Deleted";
    }
    return new Date(deletedTimestamp).toLocaleString(undefined, {
        dateStyle: "medium",
        timeStyle: "short",
    });
};

const RecycleBinDialog = ({
    open,
    onClose,
    fetchProjects,
    onDeleteAll,
    onRestore,
}: RecycleBinDialogProps) => {
    const theme = useTheme();
    const [projects, setProjects] = useState<RecycledProject[]>([]);
    const [confirmingDeleteAll, setConfirmingDeleteAll] = useState(false);
    const [projectToRestore, setProjectToRestore] = useState<RecycledProject | null>(null);

    useEffect(() => {
        if (open) {
            void fetchProjects()
                .then(setProjects)
                .catch(() => {
                    setProjects([]);
                });
        }
    }, [fetchProjects, open]);

    const handleDeleteAll = async () => {
        try {
            const result = await onDeleteAll();
            const deletedIds = new Set(result.deletedProjectIds.map(String));
            setProjects((currentProjects) =>
                currentProjects.filter((project) => !deletedIds.has(project.id)),
            );
            setConfirmingDeleteAll(false);
        } catch {
            // The hook opens the shared error dialog.
        }
    };

    const handleRestore = async () => {
        if (!projectToRestore) {
            return;
        }
        try {
            await onRestore(projectToRestore.id);
            setProjects((currentProjects) =>
                currentProjects.filter((project) => project.id !== projectToRestore.id),
            );
            setProjectToRestore(null);
        } catch {
            // The hook opens the shared error dialog.
        }
    };

    return (
        <>
            <Dialog
                open={open}
                onClose={onClose}
                fullWidth
                maxWidth="sm"
                PaperProps={{
                    sx: {
                        borderRadius: 2,
                        overflow: "hidden",
                    },
                }}
            >
                <DialogTitle sx={{ pb: 1 }}>
                    <Stack direction="row" spacing={1.5} alignItems="center">
                        <Box
                            sx={{
                                display: "flex",
                                alignItems: "center",
                                justifyContent: "center",
                                width: 40,
                                height: 40,
                                borderRadius: "50%",
                                bgcolor: alpha(theme.palette.error.main, 0.12),
                                color: "error.main",
                            }}
                        >
                            <DeleteOutlineIcon />
                        </Box>
                        <Box>
                            <Typography variant="h6" component="span">
                                Recycle bin
                            </Typography>
                            <Typography variant="body2" color="text.secondary" display="block">
                                {projects.length === 0
                                    ? "No projects waiting to be removed"
                                    : `${projects.length} project${projects.length === 1 ? "" : "s"} can be permanently deleted`}
                            </Typography>
                        </Box>
                    </Stack>
                    <DialogCloseIconButton onClose={onClose} />
                </DialogTitle>
                <Divider />
                <DialogContent sx={{ px: 3, py: 3, bgcolor: alpha(theme.palette.background.default, 0.6) }}>
                    {projects.length === 0 ? (
                        <Box
                            sx={{
                                display: "flex",
                                flexDirection: "column",
                                alignItems: "center",
                                textAlign: "center",
                                py: 4,
                                px: 2,
                            }}
                        >
                            <Box
                                sx={{
                                    display: "flex",
                                    alignItems: "center",
                                    justifyContent: "center",
                                    width: 72,
                                    height: 72,
                                    borderRadius: "50%",
                                    mb: 2,
                                    bgcolor: alpha(theme.palette.text.secondary, 0.08),
                                    color: "text.secondary",
                                }}
                            >
                                <FolderOffOutlinedIcon sx={{ fontSize: 36 }} />
                            </Box>
                            <Typography variant="subtitle1" gutterBottom>
                                Your recycle bin is empty
                            </Typography>
                            <Typography variant="body2" color="text.secondary" maxWidth={320}>
                                Projects you move to the recycle bin from the catalog will appear here
                                before they are permanently deleted.
                            </Typography>
                        </Box>
                    ) : (
                        <Stack spacing={1.5}>
                            {projects.map((project) => (
                                <Paper
                                    key={project.id}
                                    variant="outlined"
                                    sx={{
                                        px: 2,
                                        py: 1.5,
                                        borderRadius: 1.5,
                                        display: "flex",
                                        alignItems: "center",
                                        justifyContent: "space-between",
                                        gap: 2,
                                        bgcolor: "background.paper",
                                    }}
                                >
                                    <Box sx={{ minWidth: 0 }}>
                                        <Typography variant="subtitle1" noWrap>
                                            {project.name}
                                        </Typography>
                                        <Typography variant="caption" color="text.secondary">
                                            {formatDeletedAt(project.deletedTimestamp)}
                                        </Typography>
                                    </Box>
                                    {/* <Chip
                                        size="small"
                                        label="Pending deletion"
                                        color="warning"
                                        variant="outlined"
                                    /> */}
                                    <Button variant="outlined" onClick={() => setProjectToRestore(project)}>
                                        Restore
                                    </Button>
                                </Paper>
                            ))}
                        </Stack>
                    )}
                </DialogContent>
                <Divider />
                <DialogActions sx={{ px: 3, py: 2, gap: 1 }}>
                    <Button onClick={onClose} sx={{ textTransform: "none" }}>
                        Close
                    </Button>
                    <Box sx={{ flexGrow: 1 }} />
                    <Button
                        color="error"
                        variant="contained"
                        startIcon={<DeleteForeverIcon />}
                        disabled={projects.length === 0}
                        onClick={() => setConfirmingDeleteAll(true)}
                        sx={{ textTransform: "none" }}
                    >
                        Delete all permanently
                    </Button>
                </DialogActions>
            </Dialog>
        
            <Dialog
                open={confirmingDeleteAll}
                onClose={() => setConfirmingDeleteAll(false)}
                maxWidth="sm"
                fullWidth
            >
                <DialogTitle>
                    Permanently delete recycled projects?
                    <DialogCloseIconButton onClose={() => setConfirmingDeleteAll(false)} />
                </DialogTitle>
                <DialogContent dividers>
                    <Typography>
                        This permanently deletes every project in your recycle bin ({projects.length}{" "}
                        {projects.length === 1 ? "project" : "projects"}). This action is irreversible
                        and cannot be undone.
                    </Typography>
                </DialogContent>
                <DialogActions>
                    <Button onClick={() => setConfirmingDeleteAll(false)}>Cancel</Button>
                    <Button
                        color="error"
                        variant="contained"
                        startIcon={<DeleteForeverIcon />}
                        onClick={handleDeleteAll}
                    >
                        Permanently delete all
                    </Button>
                </DialogActions>
            </Dialog>
            <Dialog
                open={projectToRestore !== null}
                onClose={() => setProjectToRestore(null)}
                maxWidth="sm"
                fullWidth
            >
                <DialogTitle>
                    Restore project?
                    <DialogCloseIconButton onClose={() => setProjectToRestore(null)} />
                </DialogTitle>
                <DialogContent dividers>
                    <Typography>
                        Restore {projectToRestore ? `"${projectToRestore.name}"` : "this project"} to the dashboard?
                    </Typography>
                </DialogContent>
                <DialogActions>
                    <Button onClick={() => setProjectToRestore(null)}>Cancel</Button>
                    <Button variant="contained" onClick={handleRestore}>
                        Restore project
                    </Button>
                </DialogActions>
            </Dialog>
        </>
    );
};

export default RecycleBinDialog;
