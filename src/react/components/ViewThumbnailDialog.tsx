import {
    Box,
    Button,
    Card,
    CardActionArea,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Divider,
    IconButton,
    Typography,
    Grid2,
} from "@mui/material";
import { Close as CloseIcon } from "@mui/icons-material";
import { useCallback, useEffect, useState, type Dispatch, type SetStateAction } from "react";
import { ErrorBoundary } from "react-error-boundary";
import ErrorDisplay from "@/charts/dialogs/ErrorDisplay";

export type ViewThumbnailDialogProps = {
    open: boolean;
    setOpen: Dispatch<SetStateAction<boolean>>;
};

const ViewThumbnailDialog = ({ open, setOpen }: ViewThumbnailDialogProps) => {
    const viewManager = window.mdv.chartManager.viewManager;
    const [viewList, setViewList] = useState<{ name: string; image: string }[]>([]);

    const getViewList = useCallback(async () => {
        const list = await viewManager.getViewDetails();
        setViewList(list);
    }, [viewManager]);

    useEffect(() => {
        if (viewManager && open) {
            getViewList();
        }
    }, [viewManager, getViewList, open]);

    const onClose = () => {
        setOpen(false);
    };

    const handleCardClick = (name: string) => {
        viewManager.changeView(name);
        onClose();
    };

    return (
        <Dialog open={open} onClose={onClose} fullWidth maxWidth="md">
            <DialogTitle>
                Browse View Gallery
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
                <ErrorBoundary
                    FallbackComponent={({ error }) => (
                        <ErrorDisplay error={error} title="Error displaying view gallery" />
                    )}
                >
                    <Box sx={{ flexGrow: 1 }}>
                        <Grid2 container spacing={4}>
                            {viewList.map((view, index) => (
                                <Grid2 key={`${view.name}-${index}`} size={6}>
                                    <Card sx={{ boxShadow: 20 }}>
                                        <CardActionArea onClick={() => handleCardClick(view.name)}>
                                            <Box sx={{ display: "flex", justifyContent: "center", padding: "2%" }}>
                                                <Typography sx={{ fontWeight: "bold" }}>{view.name}</Typography>
                                            </Box>
                                            <Divider />
                                            <Box sx={{ padding: "1%" }}>
                                                <img
                                                    src={view.image}
                                                    alt={`${view.name} snapshot`}
                                                />
                                            </Box>
                                        </CardActionArea>
                                    </Card>
                                </Grid2>
                            ))}
                        </Grid2>
                    </Box>
                </ErrorBoundary>
            </DialogContent>
            <DialogActions>
                <Button onClick={onClose} color="primary" sx={{ mr: 2 }}>
                    Cancel
                </Button>
            </DialogActions>
        </Dialog>
    );
};

export default ViewThumbnailDialog;
