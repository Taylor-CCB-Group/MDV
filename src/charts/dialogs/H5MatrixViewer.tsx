import {
    Box,
    Paper,
    Stack,
    Tab,
    Table,
    TableBody,
    TableCell,
    TableContainer,
    TableHead,
    TableRow,
    Tabs,
    Typography,
} from "@mui/material";
import type React from "react";
import { useCallback, useMemo, useState } from "react";

export type Matrix = number[][];
export type MatrixRecord = Record<string, Matrix | null>;

export interface MatrixData {
    readonly X: Matrix | null;
    readonly layers: MatrixRecord;
    readonly obsm: MatrixRecord;
    readonly varm: MatrixRecord;
    readonly obsp: MatrixRecord;
    readonly varp: MatrixRecord;
}

export interface MatrixStats {
    readonly min: number;
    readonly max: number;
    readonly mean: number;
    readonly dimensions: string;
}

export interface SelectedCell {
    readonly row: number;
    readonly col: number;
    readonly value: number;
}

export interface VisibleMatrixData {
    readonly rows: ReadonlyArray<ReadonlyArray<number>>;
    readonly hasMoreRows: boolean;
    readonly hasMoreCols: boolean;
    readonly visibleColumns: ReadonlyArray<number>;
    readonly rowIndices: ReadonlyArray<number>;
    readonly colIndices: ReadonlyArray<number>;
}

export interface H5MatrixViewerProps {
    readonly matrices: Readonly<MatrixData>;
}

export type MatrixEntry = readonly [string, Matrix | null];
export type MatrixEntries = ReadonlyArray<MatrixEntry>;

const DEFAULT_VISIBLE_ITEMS = 10;

const H5MatrixViewer: React.FC<H5MatrixViewerProps> = ({ matrices }) => {
    const [activeMatrixTab, setActiveMatrixTab] = useState<number>(0);
    const [selectedCell, setSelectedCell] = useState<SelectedCell | null>(null);

    const matrixEntries = useMemo<MatrixEntries>(() => {
        // Find the first non-empty matrix category
        if (matrices.X) {
            return [["X", matrices.X]] as const;
        }
        if (Object.keys(matrices.layers).length > 0) {
            return Object.entries(matrices.layers);
        }
        if (Object.keys(matrices.obsm).length > 0) {
            return Object.entries(matrices.obsm);
        }
        if (Object.keys(matrices.varm).length > 0) {
            return Object.entries(matrices.varm);
        }
        if (Object.keys(matrices.obsp).length > 0) {
            return Object.entries(matrices.obsp);
        }
        if (Object.keys(matrices.varp).length > 0) {
            return Object.entries(matrices.varp);
        }
        return [] as MatrixEntries;
    }, [matrices]);

    const calculateStats = useCallback((matrix: Matrix): MatrixStats => {
        if (!matrix?.length) {
            return { min: 0, max: 0, mean: 0, dimensions: "0 × 0" } as const;
        }

        const chunkSize = 1000;
        let min = Number.POSITIVE_INFINITY;
        let max = Number.NEGATIVE_INFINITY;
        let sum = 0;
        let count = 0;

        for (let i = 0; i < matrix.length; i += chunkSize) {
            const chunk = matrix.slice(i, i + chunkSize);
            for (const row of chunk) {
                for (const value of row) {
                    min = Math.min(min, value);
                    max = Math.max(max, value);
                    sum += value;
                    count++;
                }
            }
        }

        return {
            min: Number(min.toFixed(4)),
            max: Number(max.toFixed(4)),
            mean: Number((sum / count).toFixed(4)),
            dimensions: `${matrix.length} × ${matrix[0].length}`,
        } as const;
    }, []);

    const currentMatrix = useMemo<MatrixEntry | null>(() => {
        return matrixEntries.length ? matrixEntries[activeMatrixTab] : null;
    }, [matrixEntries, activeMatrixTab]);

    const visibleMatrixData = useMemo((): VisibleMatrixData | null => {
        if (!currentMatrix?.[1]) return null;

        const matrix = currentMatrix[1];
        const visibleRows = Math.min(DEFAULT_VISIBLE_ITEMS, matrix.length);
        const visibleCols = Math.min(DEFAULT_VISIBLE_ITEMS, matrix[0].length);
        const rows = matrix.slice(0, visibleRows);

        return {
            rows,
            hasMoreRows: matrix.length > visibleRows,
            hasMoreCols: matrix[0].length > visibleCols,
            visibleColumns: matrix[0].slice(0, visibleCols),
            rowIndices: Array.from({ length: visibleRows }, (_, i) => i),
            colIndices: Array.from({ length: visibleCols }, (_, i) => i),
        } as const;
    }, [currentMatrix]);

    const getCellColor = useCallback(
        (value: number, min: number, max: number): string => {
            const normalizedValue = (value - min) / (max - min);
            return `rgba(0, 0, 255, ${normalizedValue * 0.2})` as const;
        },
        [],
    );

    const getEllipsisKey = useCallback(
        (
            name: string,
            type: string,
            position: string,
            matrix: Matrix,
        ): string => {
            return `${name}-${type}-ellipsis-${position}-${matrix.length}-${matrix[0].length}` as const;
        },
        [],
    );

    const renderMatrix = useCallback(
        (matrixData: MatrixEntry) => {
            const [name, matrix] = matrixData;
            if (!matrix || matrix.length === 0 || !visibleMatrixData)
                return null;

            const stats = calculateStats(matrix);
            const { rows, hasMoreRows, hasMoreCols, colIndices, rowIndices } =
                visibleMatrixData;

            return (
                <Box sx={{ width: "100%" }}>
                    <Stack
                        direction="row"
                        justifyContent="space-between"
                        alignItems="center"
                        sx={{ mb: 2 }}
                    >
                        <Stack spacing={0.5}>
                            <Typography variant="subtitle1">
                                Matrix Statistics
                            </Typography>
                            <Typography variant="body2" color="text.secondary">
                                Dimensions: {stats.dimensions} | Range: [
                                {stats.min}, {stats.max}] | Mean: {stats.mean}
                            </Typography>
                        </Stack>
                    </Stack>

                    <TableContainer
                        sx={{
                            height: "400px",
                            maxWidth: "800px",
                            overflow: "auto",
                            "& td, & th": {
                                minWidth: "100px",
                                padding: "6px 16px",
                            },
                            "& th:first-of-type, & td:first-of-type": {
                                position: "sticky",
                                left: 0,
                                backgroundColor: "background.paper",
                                zIndex: 2,
                            },
                            "& th": {
                                position: "sticky",
                                top: 0,
                                backgroundColor: "background.paper",
                                zIndex: 1,
                            },
                            "& th:first-of-type": {
                                zIndex: 3,
                            },
                        }}
                    >
                        <Table size="small" stickyHeader>
                            <TableHead>
                                <TableRow>
                                    <TableCell
                                        sx={{ bgcolor: "grey.100", width: 60 }}
                                    >
                                        #
                                    </TableCell>
                                    {colIndices.map((colIndex) => (
                                        <TableCell
                                            key={`${name}-col-${colIndex}-${matrix[0].length}`}
                                        >
                                            {colIndex}
                                        </TableCell>
                                    ))}
                                    {hasMoreCols && (
                                        <TableCell
                                            key={getEllipsisKey(
                                                name,
                                                "header",
                                                "right",
                                                matrix,
                                            )}
                                        >
                                            ...
                                        </TableCell>
                                    )}
                                </TableRow>
                            </TableHead>
                            <TableBody>
                                {rows.map((row, rowIdx) => {
                                    const actualRowIndex = rowIndices[rowIdx];
                                    return (
                                        <TableRow
                                            key={`${name}-row-${actualRowIndex}-${matrix.length}`}
                                        >
                                            <TableCell
                                                sx={{ bgcolor: "grey.100" }}
                                            >
                                                {actualRowIndex}
                                            </TableCell>
                                            {row
                                                .slice(0, DEFAULT_VISIBLE_ITEMS)
                                                .map((cell, colIdx) => {
                                                    const actualColIndex =
                                                        colIndices[colIdx];
                                                    return (
                                                        <TableCell
                                                            key={`${name}-cell-${actualRowIndex}-${actualColIndex}-${matrix.length}-${matrix[0].length}`}
                                                            sx={{
                                                                bgcolor:
                                                                    getCellColor(
                                                                        cell,
                                                                        stats.min,
                                                                        stats.max,
                                                                    ),
                                                                cursor: "pointer",
                                                                outline:
                                                                    selectedCell?.row ===
                                                                        actualRowIndex &&
                                                                    selectedCell?.col ===
                                                                        actualColIndex
                                                                        ? "2px solid #1976d2"
                                                                        : "none",
                                                                transition:
                                                                    "background-color 0.2s",
                                                            }}
                                                            onClick={() =>
                                                                setSelectedCell(
                                                                    {
                                                                        row: actualRowIndex,
                                                                        col: actualColIndex,
                                                                        value: cell,
                                                                    },
                                                                )
                                                            }
                                                        >
                                                            {cell.toFixed(4)}
                                                        </TableCell>
                                                    );
                                                })}
                                            {hasMoreCols && (
                                                <TableCell
                                                    key={getEllipsisKey(
                                                        name,
                                                        "row",
                                                        `${actualRowIndex}-right`,
                                                        matrix,
                                                    )}
                                                >
                                                    ...
                                                </TableCell>
                                            )}
                                        </TableRow>
                                    );
                                })}
                                {hasMoreRows && (
                                    <TableRow
                                        key={getEllipsisKey(
                                            name,
                                            "footer",
                                            "row",
                                            matrix,
                                        )}
                                    >
                                        <TableCell
                                            key={getEllipsisKey(
                                                name,
                                                "footer",
                                                "index",
                                                matrix,
                                            )}
                                        >
                                            ...
                                        </TableCell>
                                        {Array.from({
                                            length: hasMoreCols
                                                ? DEFAULT_VISIBLE_ITEMS + 1
                                                : DEFAULT_VISIBLE_ITEMS,
                                        }).map((_, colIdx) => (
                                            <TableCell
                                                key={getEllipsisKey(
                                                    name,
                                                    "footer",
                                                    `cell-${matrix.length}-${colIdx}`,
                                                    matrix,
                                                )}
                                            >
                                                ...
                                            </TableCell>
                                        ))}
                                    </TableRow>
                                )}
                            </TableBody>
                        </Table>
                    </TableContainer>
                </Box>
            );
        },
        [
            visibleMatrixData,
            calculateStats,
            getCellColor,
            selectedCell,
            getEllipsisKey,
        ],
    );

    return (
        <Box sx={{ width: "100%", mt: 2 }}>
            <Paper elevation={4}>
                <Box sx={{ p: 2 }}>
                    {matrixEntries.length > 0 && (
                        <>
                            <Tabs
                                value={activeMatrixTab}
                                onChange={(_, newValue) =>
                                    setActiveMatrixTab(newValue)
                                }
                                variant="scrollable"
                                scrollButtons="auto"
                                sx={{
                                    mb: 2,
                                    borderBottom: 1,
                                    borderColor: "divider",
                                }}
                            >
                                {matrixEntries.map(([name]) => (
                                    <Tab key={name} label={name} />
                                ))}
                            </Tabs>
                            {currentMatrix && renderMatrix(currentMatrix)}
                        </>
                    )}
                </Box>
            </Paper>
        </Box>
    );
};

export default H5MatrixViewer;
