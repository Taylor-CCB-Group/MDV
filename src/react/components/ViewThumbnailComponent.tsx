import GridViewIcon from "@mui/icons-material/GridView";
import { observer } from "mobx-react-lite";
import { useEffect, useState } from "react";
import ViewThumbnailDialog from "./ViewThumbnailDialog";
import IconWithTooltip from "./IconWithTooltip";

const ViewThumbnailComponent = observer(() => {
    const viewManager = window.mdv.chartManager.viewManager;
    const [open, setOpen] = useState(false);

    // Auto-open gallery when viewManager.showGallery is set (e.g. on project load with gallery default)
    useEffect(() => {
        if (viewManager.showGallery && !open) {
            setOpen(true);
            viewManager.setShowGallery(false);
        }
    }, [viewManager.showGallery, open, viewManager]);

    const onOpen = () => {
        setOpen(true);
    };

    return (
        <>
            <IconWithTooltip tooltipText="Browse View Gallery" onClick={onOpen}>
                <GridViewIcon />
            </IconWithTooltip>

            <ViewThumbnailDialog open={open} setOpen={setOpen} />
        </>
    );
});

export default ViewThumbnailComponent;
