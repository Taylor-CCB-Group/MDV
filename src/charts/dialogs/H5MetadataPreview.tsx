import { Box, Grid, Paper, Typography } from "@mui/material";
import React, { useState, useEffect } from "react";
import H5JsonViewer from "./H5JsonViewer";
import H5MatrixViewer from "./H5MatrixViewer";

export interface H5MetadataPreviewProps {
    metadata: {
        uns: Record<string, any>;
        obs: Record<string, any>;
        var: Record<string, any>;
        X: number[][] | null;
        layers: Record<string, number[][] | null>;
        obsm: Record<string, number[][] | null>;
        varm: Record<string, number[][] | null>;
        obsp: Record<string, number[][] | null>;
        varp: Record<string, number[][] | null>;
    };
}

const H5MetadataPreview = ({ metadata }: H5MetadataPreviewProps) => {
    const [selectedMetadata, setSelectedMetadata] = useState("");
    const [availableMetadata, setAvailableMetadata] = useState<
        Array<{ key: string; label: string; hasData: boolean }>
    >([]);

    useEffect(() => {
        const available = [
            {
                key: "uns",
                label: "Unstructured Data (uns)",
                hasData: !!metadata.uns && Object.keys(metadata.uns).length > 0,
            },
            {
                key: "obs",
                label: "Observation Data (obs)",
                hasData: !!metadata.obs && Object.keys(metadata.obs).length > 0,
            },
            {
                key: "var",
                label: "Variable Data (var)",
                hasData: !!metadata.var && Object.keys(metadata.var).length > 0,
            },
            { key: "X", label: "Expression Matrix (X)", hasData: !!metadata.X },
            {
                key: "layers",
                label: "Layers",
                hasData:
                    !!metadata.layers &&
                    Object.keys(metadata.layers).length > 0,
            },
            {
                key: "obsm",
                label: "Observation Matrices (obsm)",
                hasData:
                    !!metadata.obsm && Object.keys(metadata.obsm).length > 0,
            },
            {
                key: "varm",
                label: "Variable Matrices (varm)",
                hasData:
                    !!metadata.varm && Object.keys(metadata.varm).length > 0,
            },
            {
                key: "obsp",
                label: "Observation Pairs (obsp)",
                hasData:
                    !!metadata.obsp && Object.keys(metadata.obsp).length > 0,
            },
            {
                key: "varp",
                label: "Variable Pairs (varp)",
                hasData:
                    !!metadata.varp && Object.keys(metadata.varp).length > 0,
            },
        ];

        setAvailableMetadata(available);

        const firstAvailable = available.find((m) => m.hasData);
        if (firstAvailable) {
            setSelectedMetadata(firstAvailable.key);
        }
    }, [metadata]);

    const isMatrixData = (key: string) =>
        ["X", "layers", "obsm", "varm", "obsp", "varp"].includes(key);

    return (
        <Box
            sx={{
                width: "865px",
                height: "710px",
                px: 2,
            }}
        >
            <Paper
                elevation={4}
                sx={{
                    overflowX: "auto",
                    border: "1px solid",
                    borderColor: "var(--main_panel_color)",
                    position: "relative", // Added to establish stacking context
                }}
            >
                <Paper
                    elevation={1}
                    sx={{
                        position: "sticky",
                        top: 0,
                        zIndex: 1,
                        borderColor: "var(--main_panel_color)",
                        bgcolor: "var(--background_color)",
                    }}
                >
                    <Typography
                        variant="h6"
                        sx={{
                            p: 2,
                            textAlign: "center",
                            backgroundColor: "var(--background_color)",
                        }}
                    >
                        H5 File Metadata
                    </Typography>

                    <Box sx={{ px: 2, pb: 2 }}>
                        <select
                            value={selectedMetadata}
                            onChange={(e) =>
                                setSelectedMetadata(e.target.value)
                            }
                            className="block w-full px-3 py-2 border border-gray-300 bg-white rounded-md shadow-sm focus:outline-none focus:ring-blue-500 focus:border-blue-500 sm:text-sm dark:bg-gray-800 dark:border-gray-600 dark:text-white"
                        >
                            {availableMetadata.map(
                                ({ key, label, hasData }) => (
                                    <option
                                        key={key}
                                        value={key}
                                        disabled={!hasData}
                                        className={
                                            !hasData ? "text-gray-400" : ""
                                        }
                                    >
                                        {label} {!hasData ? "(empty)" : ""}
                                    </option>
                                ),
                            )}
                        </select>
                    </Box>
                </Paper>

                <Box
                    sx={{
                        p: 2,
                        backgroundColor: "var(--primary_background)",
                        position: "relative", // For proper stacking
                        zIndex: 0, // Below the header
                        height: isMatrixData(selectedMetadata)
                            ? "auto"
                            : "600px",
                        display: "flex",
                        flexDirection: "column",
                    }}
                >
                    {selectedMetadata && !isMatrixData(selectedMetadata) && (
                        <Paper
                            elevation={0}
                            sx={{
                                p: 2,
                                backgroundColor: "white",
                                flexGrow: 1,
                                overflowY: "auto",
                                overflowX: "auto",
                                border: "1px solid",
                                borderColor: "rgb(209, 213, 219)",
                                height: "500px",
                            }}
                        >
                            <H5JsonViewer
                                data={
                                    metadata[
                                        selectedMetadata as keyof typeof metadata
                                    ] || {}
                                }
                                title={
                                    availableMetadata.find(
                                        (m) => m.key === selectedMetadata,
                                    )?.label || ""
                                }
                                initiallyExpanded={true}
                                maxHeight="none"
                            />
                        </Paper>
                    )}

                    {selectedMetadata && isMatrixData(selectedMetadata) && (
                        <Box>
                            <H5MatrixViewer
                                matrices={{
                                    X:
                                        selectedMetadata === "X"
                                            ? metadata.X
                                            : null,
                                    layers:
                                        selectedMetadata === "layers"
                                            ? metadata.layers
                                            : {},
                                    obsm:
                                        selectedMetadata === "obsm"
                                            ? metadata.obsm
                                            : {},
                                    varm:
                                        selectedMetadata === "varm"
                                            ? metadata.varm
                                            : {},
                                    obsp:
                                        selectedMetadata === "obsp"
                                            ? metadata.obsp
                                            : {},
                                    varp:
                                        selectedMetadata === "varp"
                                            ? metadata.varp
                                            : {},
                                }}
                            />
                        </Box>
                    )}
                </Box>
            </Paper>
        </Box>
    );
};

export default H5MetadataPreview;
