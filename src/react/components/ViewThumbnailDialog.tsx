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
import { Clear as ClearIcon, Close as CloseIcon, Image as ImageIcon } from "@mui/icons-material";
import { useCallback, useEffect, useState, type Dispatch, type SetStateAction } from "react";
import { ErrorBoundary } from "react-error-boundary";
import DebugErrorComponent from "@/charts/dialogs/DebugErrorComponent";

export type ViewThumbnailDialogProps = {
    open: boolean;
    setOpen: Dispatch<SetStateAction<boolean>>;
};

export type ViewImageComponentProps = {
    imgSrc: string;
    viewName: string;
};

export const ViewImageComponent = ({imgSrc, viewName}: ViewImageComponentProps) => {
    const [hasError, setHasError] = useState(!imgSrc);

    return !hasError ? (
        <img
            src={imgSrc}
            alt={`${viewName} snapshot`}
            style={{ objectFit: 'contain', aspectRatio: 3/2 }}
            onError={() => setHasError(true)}
        />
    ) : (
        <Box sx={{ 
            width: "100%", 
            height: "100%", 
            aspectRatio: 3/2,
            display: "flex",
            justifyContent: "center",
            alignItems: "center", 
        }}>
            <ImageIcon 
                sx={{ 
                    fontSize: "5rem",
                    color: "text.secondary"
                }} 
            />
        </Box>
    );
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

    const handleCardClick = async (name: string) => {
        viewManager.checkAndChangeView(name);
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
            <DialogContent dividers sx={{height: "95vh"}}>
                <ErrorBoundary
                    FallbackComponent={({ error }) => (
                        <DebugErrorComponent error={error} title="Error displaying view gallery" />
                    )}
                >
                    {loading ? (
                        <Box
                            sx={{
                                display: "flex",
                                justifyContent: "center",
                                alignItems: "center",
                                height: "100%",
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
                                                <Card sx={{ boxShadow: 20, bgcolor: "inherit" }}>
                                                    <CardActionArea onClick={() => handleCardClick(view.name)}>
                                                        <Box sx={{ padding: "1%", bgcolor: "background.default" }}>
                                                            <ViewImageComponent imgSrc={view.image} viewName={view.name} />
                                                        </Box>
                                                        <Divider />
                                                        <Box
                                                            sx={{
                                                                display: "flex",
                                                                justifyContent: "center",
                                                                paddingY: "3%",
                                                            }}
                                                        >
                                                            <Typography sx={{ fontWeight: "bold" }}>
                                                                {view.name}
                                                            </Typography>
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
