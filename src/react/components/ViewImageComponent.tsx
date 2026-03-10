import { Box } from "@mui/material";
import { useEffect, useState } from "react";
import ImageIcon from "@mui/icons-material/Image";

export type ViewImageComponentProps = {
    imgSrc: string;
    viewName: string;
};

const ViewImageComponent = ({ imgSrc, viewName }: ViewImageComponentProps) => {
    const [hasError, setHasError] = useState(!imgSrc);

    useEffect(() => {
        // update hasError when imgSrc updates
        setHasError(!imgSrc);
    }, [imgSrc]);

    const ZOOM_LEVEL = 1.75;
    const VERTICAL_OFFSET_PERCENT = 25;
    const baseWidth = `${100 / ZOOM_LEVEL}%`;

    return (
        <Box
            sx={{
                width: "100%",
                height: "100%",
                aspectRatio: 3 / 2,
                display: "flex",
                justifyContent: "center",
                alignItems: "flex-start",
                overflow: "hidden",
            }}
        >
            {!hasError ? (
                <img
                    src={imgSrc}
                    alt={`${viewName} snapshot`}
                    style={{
                        width: baseWidth,
                        height: "auto",
                        transform: `scale(${ZOOM_LEVEL}) translateY(${VERTICAL_OFFSET_PERCENT}%)`,
                        transition: "transform 0.2s ease-in-out",
                        transformOrigin: "center top",
                    }}
                    onError={() => setHasError(true)}
                />
            ) : (
                <ImageIcon
                    sx={{
                        fontSize: "5rem",
                        color: "text.secondary",
                    }}
                />
            )}
        </Box>
    );
};

export default ViewImageComponent;