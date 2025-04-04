import {
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
} from "@mui/material";
import { Close as CloseIcon } from "@mui/icons-material";
import { useCallback, useEffect, useState, type Dispatch, type SetStateAction } from "react";

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
        if (viewManager) {
            getViewList();
        } 
    }, [viewManager, getViewList]);

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
                <div style={{ display: "flex", flexWrap: "wrap", gap: 16 }}>
                    {viewList.map((view, index) => (
                        <Card key={`${view.name}-${index}`}>
                            <CardActionArea onClick={() => handleCardClick(view.name)}>
                                <img src={view.image} alt="" />
                                <Divider />
                                <div style={{ display: "flex", justifyContent: "center", padding: "3%" }}>
                                    <Typography color="text.secondary">{view.name}</Typography>
                                </div>
                            </CardActionArea>
                        </Card>
                    ))}
                </div>
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
