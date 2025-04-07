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
    CircularProgress,
    TextField,
    InputAdornment,
} from "@mui/material";
import { Clear as ClearIcon, Close as CloseIcon } from "@mui/icons-material";
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
    const [filteredViewList, setFilteredViewList] = useState<{ name: string; image: string }[]>([]);
    const [loading, setLoading] = useState(false);
    const [viewName, setViewName] = useState<string>();

    const getViewList = useCallback(async () => {
        setLoading(true);
        const list = await viewManager.getViewDetails();
        setViewList(list);
        setFilteredViewList(list);
        setLoading(false);
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

    const onInputChange = (inputText: string) => {
        setViewName(inputText);
        const input = inputText.toLowerCase().split(" ");
        const tempList = viewList.filter((view) => {
            const name = view.name.toLowerCase();
            // if any of the input words are not in the name, return false
            if (!input.some(i => !name.includes(i))) return view;
        });
        setFilteredViewList(tempList);
    };

    const handleClearInput = () => {
        setViewName("");
        setFilteredViewList(viewList);
    };

    return (
        <Dialog open={open} onClose={onClose} fullWidth maxWidth="xl">
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
                    {loading ? (
                        <Box
                            sx={{
                                display: "flex",
                                justifyContent: "center",
                                alignItems: "center",
                                height: "20vh",
                                width: "100%",
                            }}
                        >
                            <CircularProgress />
                        </Box>
                    ) : (
                        <Box sx={{ flexGrow: 1 }}>
                            <Box sx={{ mt: 1, mb: 3 }}>
                                <TextField
                                    fullWidth
                                    placeholder="Enter the view name to filter the views"
                                    value={viewName}
                                    onChange={(e) => onInputChange(e.target.value)}
                                    slotProps={{
                                        input: {
                                            endAdornment: (
                                                <InputAdornment position="end">
                                                    <IconButton onClick={handleClearInput}>
                                                        <ClearIcon />
                                                    </IconButton>
                                                </InputAdornment>
                                            )
                                        }
                                    }}
                                />
                            </Box>
                            {filteredViewList.length === 0 ? (
                                <Box
                                    sx={{
                                        display: "flex",
                                        justifyContent: "center",
                                        alignItems: "center",
                                        height: "20vh",
                                        width: "100%",
                                    }}
                                >
                                    <Typography variant="h6">No views found</Typography>
                                </Box>
                            ) : (
                                <>
                                    <Grid2 container spacing={4}>
                                        {filteredViewList.map((view, index) => (
                                            <Grid2 key={`${view.name}-${index}`} size={{md: 4, sm: 6}}>
                                                <Card sx={{ boxShadow: 20 }}>
                                                    <CardActionArea onClick={() => handleCardClick(view.name)}>
                                                        <Box
                                                            sx={{
                                                                display: "flex",
                                                                justifyContent: "center",
                                                                padding: "2%",
                                                            }}
                                                        >
                                                            <Typography sx={{ fontWeight: "bold" }}>
                                                                {view.name}
                                                            </Typography>
                                                        </Box>
                                                        <Divider />
                                                        <Box sx={{ padding: "1%" }}>
                                                            <img
                                                                src={view.image}
                                                                alt={`${view.name} snapshot`}
                                                                style={{ objectFit: 'inherit', aspectRatio: 4/3 }}
                                                            />
                                                        </Box>
                                                    </CardActionArea>
                                                </Card>
                                            </Grid2>
                                        ))}
                                    </Grid2>
                                </>
                            )}
                        </Box>
                    )}
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
