import GridViewIcon from "@mui/icons-material/GridView";
import { observer } from "mobx-react-lite";
import { useEffect, useState } from "react";
import ViewThumbnailDialog from "./ViewThumbnailDialog";
import IconWithTooltip from "./IconWithTooltip";
import { useChartManager } from "../hooks";

const ViewThumbnailComponent = observer(() => {
    const cm = useChartManager();
    const viewManager = cm.viewManager;
    const [open, setOpen] = useState(false);

    // Auto-open gallery when viewManager.showGallery is set (e.g. on project load with gallery default)
    useEffect(() => {
        if (!viewManager.showGallery) return;

        if (!open) {
            setOpen(true);
        }
        // Always consume so it doesn't reopen later
        viewManager.setShowGallery(false);
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
