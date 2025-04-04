import { BrowseGallery as BrowseGalleryIcon } from "@mui/icons-material";
import { useState } from "react";
import ViewThumbnailDialog from "./ViewThumbnailDialog";
import IconWithTooltip from "./IconWithTooltip";

const ViewThumbnailComponent = () => {
    const [open, setOpen] = useState(false);

    const onOpen = () => {
        setOpen(true);
    };

    return (
        <>
            <IconWithTooltip tooltipText="Browse View Gallery" onClick={onOpen}>
                <BrowseGalleryIcon />
            </IconWithTooltip>

            <ViewThumbnailDialog open={open} setOpen={setOpen} />
        </>
    );
};

export default ViewThumbnailComponent;
