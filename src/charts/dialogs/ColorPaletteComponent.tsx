import SearchRoundedIcon from "@mui/icons-material/SearchRounded";
import {
    Box,
    Button,
    Chip,
    Divider,
    InputAdornment,
    MenuItem,
    Paper,
    Stack,
    TextField,
    Tooltip,
    Typography,
} from "@mui/material";
import { type ChangeEvent, useEffect, useMemo, useState } from "react";
import CustomPaletteDialog from "./CustomPaletteDialog";

interface ColumnOption {
    field: string;
    name: string;
}

interface ColorEntry {
    color: string;
    index: number;
    label: string;
}

interface DataSourceLike {
    getColumnColors: (column: string) => string[];
    getColumnList: (datatype: string) => ColumnOption[];
    getColumnValues: (column: string) => Array<string | number | null | undefined>;
    name: string;
}

interface ColorPaletteComponentProps {
    dataSource: DataSourceLike;
    dataSourceLabel: string;
    onApply: (column: string, colors: string[]) => void;
}

function normalizeHexColor(value: string): string | null {
    const compact = value.trim().replaceAll('"', "");
    if (!compact) {
        return null;
    }

    const match = compact.match(/^#?([0-9a-f]{3}|[0-9a-f]{6})$/i);
    if (!match) {
        return null;
    }

    const hex = match[1].toUpperCase();
    if (hex.length === 3) {
        return `#${hex
            .split("")
            .map((char) => `${char}${char}`)
            .join("")}`;
    }
    return `#${hex}`;
}

function parseHexColors(text: string) {
    const tokens = text
        .replaceAll('"', " ")
        .split(/[,\s]+/)
        .map((token) => token.trim())
        .filter(Boolean);

    const colors: string[] = [];
    let ignored = 0;

    for (const token of tokens) {
        const normalized = normalizeHexColor(token);
        if (normalized) {
            colors.push(normalized);
        } else {
            ignored += 1;
        }
    }

    return { colors, ignored, total: tokens.length };
}

function buildEntries(dataSource: DataSourceLike, column: string): ColorEntry[] {
    const values = dataSource.getColumnValues(column);
    const colors = dataSource.getColumnColors(column);

    return values
        .map((value, index) => ({
            color: normalizeHexColor(colors[index] ?? "") ?? "#808080",
            index,
            label: String(value ?? ""),
        }))
        .sort((left, right) =>
            left.label.localeCompare(right.label, undefined, {
                numeric: true,
                sensitivity: "base",
            }),
        );
}

function cloneEntries(entries: ColorEntry[]) {
    return entries.map((entry) => ({ ...entry }));
}

const scrollAreaSx = {
    flex: 1,
    minHeight: 0,
    overflowY: "auto",
    overscrollBehavior: "contain",
    scrollbarGutter: "stable",
    scrollbarWidth: "thin",
    WebkitOverflowScrolling: "touch",
    "&::-webkit-scrollbar": {
        width: 12,
    },
    "&::-webkit-scrollbar-thumb": {
        backgroundColor: "rgba(120, 120, 120, 0.45)",
        border: "2px solid transparent",
        borderRadius: 999,
        backgroundClip: "padding-box",
    },
    "&::-webkit-scrollbar-track": {
        backgroundColor: "transparent",
    },
} as const;

function ColorPaletteComponent({ dataSource, dataSourceLabel, onApply }: ColorPaletteComponentProps) {
    const columns = useMemo(() => dataSource.getColumnList("string"), [dataSource]);
    const [selectedField, setSelectedField] = useState(columns[0]?.field ?? "");
    const [entries, setEntries] = useState<ColorEntry[]>([]);
    const [initialEntries, setInitialEntries] = useState<ColorEntry[]>([]);
    const [categoryFilter, setCategoryFilter] = useState("");
    const [schemeInput, setSchemeInput] = useState("");
    const [isPaletteDialogOpen, setIsPaletteDialogOpen] = useState(false);

    useEffect(() => {
        if (columns.length === 0) {
            if (selectedField) {
                setSelectedField("");
            }
            return;
        }

        const hasSelectedField = columns.some((column) => column.field === selectedField);
        if (!hasSelectedField) {
            setSelectedField(columns[0].field);
        }
    }, [columns, selectedField]);

    useEffect(() => {
        if (!selectedField) {
            setEntries([]);
            setInitialEntries([]);
            return;
        }

        const nextEntries = buildEntries(dataSource, selectedField);
        setEntries(nextEntries);
        setInitialEntries(cloneEntries(nextEntries));
        setCategoryFilter("");
    }, [dataSource, selectedField]);

    const parsedScheme = useMemo(() => parseHexColors(schemeInput), [schemeInput]);
    const filteredEntries = useMemo(() => {
        const needle = categoryFilter.trim().toLocaleLowerCase();
        if (!needle) {
            return entries;
        }
        return entries.filter((entry) => entry.label.toLocaleLowerCase().includes(needle));
    }, [categoryFilter, entries]);

    const uniquePalette = useMemo(() => {
        const seen = new Set<string>();
        const palette: string[] = [];
        for (const entry of entries) {
            const normalized = entry.color.toUpperCase();
            if (!seen.has(normalized)) {
                seen.add(normalized);
                palette.push(normalized);
            }
        }
        return palette;
    }, [entries]);

    const isDirty = useMemo(() => {
        if (entries.length !== initialEntries.length) {
            return true;
        }
        return entries.some((entry, index) => entry.color !== initialEntries[index]?.color);
    }, [entries, initialEntries]);

    const selectedColumnName = columns.find((column) => column.field === selectedField)?.name ?? selectedField;

    const handleEntryColorChange = (entryIndex: number, color: string) => {
        setEntries((currentEntries) =>
            currentEntries.map((entry) => (entry.index === entryIndex ? { ...entry, color } : entry)),
        );
    };

    const handleApplyScheme = () => {
        if (parsedScheme.colors.length === 0) {
            return;
        }

        setEntries((currentEntries) =>
            currentEntries.map((entry, index) => ({
                ...entry,
                color: parsedScheme.colors[index % parsedScheme.colors.length],
            })),
        );
        setIsPaletteDialogOpen(false);
    };

    const handleResetField = () => {
        setEntries(cloneEntries(initialEntries));
    };

    const handleLoadCurrentPalette = () => {
        setSchemeInput(uniquePalette.join("\n"));
    };

    const handleApply = () => {
        if (!selectedField || entries.length === 0) {
            return;
        }

        const colors = new Array(entries.length);
        for (const entry of entries) {
            colors[entry.index] = entry.color;
        }

        onApply(selectedField, colors);
        setInitialEntries(cloneEntries(entries));
    };

    return (
        <>
            <Box
                sx={{
                    display: "flex",
                    flex: 1,
                    flexDirection: "column",
                    gap: 2,
                    height: "100%",
                    minHeight: 0,
                    overflow: "hidden",
                    p: { xs: 1, sm: 2 },
                }}
            >
                <Paper
                    elevation={0}
                    sx={{
                        border: 1,
                        borderColor: "divider",
                        borderRadius: 3,
                        flexShrink: 0,
                        p: { xs: 1.25, sm: 2 },
                        backgroundColor: "action.hover",
                    }}
                >
                    {/* <Typography variant="h6">Manage Color Palette</Typography> */}
                    <Typography color="text.secondary" variant="body2">
                        Choose a column and adjust colors in '{dataSourceLabel}'
                    </Typography>
                </Paper>

                {columns.length === 0 ? (
                    <Paper
                        elevation={0}
                        sx={{
                            alignItems: "center",
                            border: 1,
                            borderColor: "divider",
                            borderRadius: 3,
                            display: "flex",
                            flex: 1,
                            justifyContent: "center",
                            minHeight: 0,
                            p: 3,
                            textAlign: "center",
                        }}
                    >
                        <Box>
                            <Typography variant="subtitle1">No categorical columns are available</Typography>
                            <Typography color="text.secondary" sx={{ mt: 1 }} variant="body2">
                                This dialog works with string columns that already expose category values.
                            </Typography>
                        </Box>
                    </Paper>
                ) : (
                    <Paper
                        elevation={0}
                        sx={{
                            border: 1,
                            borderColor: "divider",
                            borderRadius: 3,
                            display: "flex",
                            flex: 1,
                            flexDirection: "column",
                            minHeight: 0,
                            overflow: "hidden",
                        }}
                    >
                        <Stack
                            spacing={2}
                            sx={{
                                flexShrink: 0,
                                p: { xs: 1.25, sm: 2.5 },
                            }}
                        >
                            <TextField
                                select
                                disabled={columns.length === 0}
                                fullWidth
                                label="Column"
                                onChange={(event) => {
                                    setSelectedField(event.target.value);
                                }}
                                value={selectedField}
                            >
                                {columns.map((column) => (
                                    <MenuItem key={column.field} value={column.field}>
                                        {column.name}
                                    </MenuItem>
                                ))}
                            </TextField>

                            <TextField
                                fullWidth
                                label="Search"
                                onChange={(event) => {
                                    setCategoryFilter(event.target.value);
                                }}
                                placeholder="Type to filter categories"
                                value={categoryFilter}
                                slotProps={{
                                    input: {
                                        startAdornment: (
                                            <InputAdornment position="start">
                                                <SearchRoundedIcon fontSize="small" />
                                            </InputAdornment>
                                        ),
                                    },
                                }}
                            />

                            <Stack
                                direction={{ sm: "row", xs: "column" }}
                                spacing={1}
                                sx={{ alignItems: { sm: "center", xs: "flex-start" }, justifyContent: "space-between" }}
                            >
                                <Stack direction="row" spacing={1} useFlexGap sx={{ flexWrap: "wrap" }}>
                                    <Chip label={`${entries.length} categories`} size="small" />
                                    <Chip
                                        label={
                                            filteredEntries.length === entries.length
                                                ? `${entries.length} shown`
                                                : `${filteredEntries.length} of ${entries.length} shown`
                                        }
                                        size="small"
                                        variant="outlined"
                                    />
                                </Stack>
                            </Stack>
                        </Stack>

                        <Divider />

                        <Box
                            sx={{
                                ...scrollAreaSx,
                                display: "flex",
                                flexDirection: "column",
                                gap: 1,
                                minHeight: 120,
                                p: { xs: 1, sm: 2 },
                            }}
                        >
                            {filteredEntries.length === 0 ? (
                                <Box
                                    sx={{
                                        alignItems: "center",
                                        color: "text.secondary",
                                        display: "flex",
                                        flex: 1,
                                        justifyContent: "center",
                                        minHeight: 220,
                                        textAlign: "center",
                                    }}
                                >
                                    <Typography variant="body2">No categories match the current search.</Typography>
                                </Box>
                            ) : (
                                filteredEntries.map((entry) => (
                                    <Box
                                        key={`${entry.label}-${entry.index}`}
                                        sx={{
                                            alignItems: "center",
                                            backgroundColor: "action.hover",
                                            borderRadius: 2,
                                            display: "grid",
                                            gap: 1.25,
                                            gridTemplateColumns: {
                                                xs: "minmax(0, 1fr) auto",
                                                sm: "minmax(0, 1fr) auto auto",
                                            },
                                            p: 1.25,
                                        }}
                                    >
                                        <Stack
                                            direction="row"
                                            spacing={1.25}
                                            sx={{ alignItems: "center", minWidth: 0 }}
                                        >
                                            <Box
                                                sx={{
                                                    backgroundColor: entry.color,
                                                    border: 1,
                                                    borderColor: "divider",
                                                    borderRadius: "50%",
                                                    flexShrink: 0,
                                                    height: 12,
                                                    width: 12,
                                                }}
                                            />
                                            <Typography noWrap variant="body2">
                                                {entry.label || "(empty value)"}
                                            </Typography>
                                        </Stack>

                                        <Chip
                                            label={entry.color.toUpperCase()}
                                            size="small"
                                            sx={{
                                                display: { xs: "none", sm: "inline-flex" },
                                                fontFamily:
                                                    'ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace',
                                            }}
                                            variant="outlined"
                                        />

                                        <Box
                                            component="input"
                                            onChange={(event: ChangeEvent<HTMLInputElement>) => {
                                                handleEntryColorChange(entry.index, event.target.value.toUpperCase());
                                            }}
                                            sx={{
                                                backgroundColor: "transparent",
                                                border: "none",
                                                cursor: "pointer",
                                                height: 36,
                                                justifySelf: "end",
                                                p: 0,
                                                width: 44,
                                                "&::-webkit-color-swatch-wrapper": {
                                                    p: 0,
                                                },
                                                "&::-webkit-color-swatch": {
                                                    border: "none",
                                                    borderRadius: "10px",
                                                },
                                                "&::-moz-color-swatch": {
                                                    border: "none",
                                                    borderRadius: "10px",
                                                },
                                            }}
                                            type="color"
                                            value={entry.color}
                                        />
                                    </Box>
                                ))
                            )}
                        </Box>

                        <Divider />

                        <Stack
                            direction={{ sm: "row", xs: "column" }}
                            spacing={1.5}
                            sx={{
                                alignItems: { sm: "center", xs: "stretch" },
                                backgroundColor: "background.paper",
                                bottom: 0,
                                flexShrink: 0,
                                justifyContent: "space-between",
                                p: { xs: 1, sm: 2 },
                                position: "sticky",
                                zIndex: 1,
                            }}
                        >
                            <Stack
                                direction={{ sm: "row", xs: "column" }}
                                spacing={1}
                                sx={{ alignItems: { sm: "center", xs: "stretch" } }}
                            >
                                <Button
                                    onClick={() => {
                                        setIsPaletteDialogOpen(true);
                                    }}
                                    variant="outlined"
                                >
                                    Custom palette
                                </Button>
                                <Button disabled={!isDirty} onClick={handleResetField} variant="text" color="error">
                                    Reset
                                </Button>
                            </Stack>

                            <Stack
                                direction="row"
                                spacing={1}
                                sx={{
                                    justifyContent: "flex-end",
                                    ml: { sm: "auto", xs: 0 },
                                    width: { xs: "100%", sm: "auto" },
                                }}
                            >
                                <Button
                                    disabled={!selectedField || entries.length === 0}
                                    onClick={handleApply}
                                    variant="contained"
                                >
                                    Apply colors
                                </Button>
                            </Stack>
                        </Stack>
                    </Paper>
                )}
            </Box>

            <CustomPaletteDialog
                open={isPaletteDialogOpen}
                parsedScheme={parsedScheme}
                schemeInput={schemeInput}
                selectedColumnName={selectedColumnName}
                onApply={handleApplyScheme}
                onClose={() => {
                    setSchemeInput("");
                    setIsPaletteDialogOpen(false);
                }}
                onLoadCurrentPalette={handleLoadCurrentPalette}
                onSchemeInputChange={setSchemeInput}
            />
        </>
    );
}

export default ColorPaletteComponent;
