import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";
import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    Paper,
    Stack,
    TextField,
    Tooltip,
    Typography,
} from "@mui/material";

interface CustomPaletteDialogProps {
    open: boolean;
    schemeInput: string;
    selectedColumnName: string;
    parsedScheme: {
        colors: string[];
    };
    onApply: () => void;
    onClose: () => void;
    onLoadCurrentPalette: () => void;
    onSchemeInputChange: (value: string) => void;
}

function buildPreviewColors(colors: string[]) {
    const colorCounts = new Map<string, number>();

    return colors.slice(0, 20).map((color) => {
        const duplicateCount = colorCounts.get(color) ?? 0;
        colorCounts.set(color, duplicateCount + 1);

        return {
            color,
            key: duplicateCount === 0 ? color : `${color}-${duplicateCount}`,
        };
    });
}

function CustomPaletteDialog({
    open,
    schemeInput,
    selectedColumnName,
    parsedScheme,
    onApply,
    onClose,
    onLoadCurrentPalette,
    onSchemeInputChange,
}: CustomPaletteDialogProps) {
    const previewColors = buildPreviewColors(parsedScheme.colors);

    return (
        <Dialog fullWidth maxWidth="sm" open={open} onClose={onClose}>
            <DialogTitle>
                Custom palette
                <DialogCloseIconButton onClose={onClose} />
            </DialogTitle>
            <DialogContent dividers>
                <Stack spacing={2}>
                    <Typography color="text.secondary" variant="body2">
                        Paste a list of hex colors. The palette will be applied across all categories in column{" "}
                        <b>{selectedColumnName}</b>.
                    </Typography>

                    <TextField
                        fullWidth
                        minRows={10}
                        multiline
                        onChange={(event) => {
                            onSchemeInputChange(event.target.value);
                        }}
                        placeholder={"#0EA5E9\n#14B8A6\n#F59E0B\n#EF4444"}
                        value={schemeInput}
                    />

                    {parsedScheme.colors.length > 0 ? (
                        <Paper
                            elevation={0}
                            sx={{
                                backgroundColor: "action.hover",
                                border: 1,
                                borderColor: "divider",
                                borderRadius: 2.5,
                                p: 1.5,
                            }}
                        >
                            <Stack direction="row" spacing={0.75} useFlexGap sx={{ flexWrap: "wrap" }}>
                                {previewColors.map(({ color, key }) => (
                                    <Tooltip key={key} title={color}>
                                        <Box
                                            sx={{
                                                backgroundColor: color,
                                                border: 1,
                                                borderColor: "divider",
                                                borderRadius: 1.5,
                                                height: 22,
                                                width: 22,
                                            }}
                                        />
                                    </Tooltip>
                                ))}
                            </Stack>
                        </Paper>
                    ) : null}
                </Stack>
            </DialogContent>
            <DialogActions
                sx={{
                    justifyContent: "space-between",
                    flexWrap: "wrap",
                    gap: 1,
                    px: 3,
                    py: 2,
                }}
            >
                <Button onClick={onLoadCurrentPalette} variant="outlined">
                    Load current palette
                </Button>
                <Stack direction="row" spacing={1}>
                    <Button onClick={onClose} variant="text" color="error">
                        Cancel
                    </Button>
                    <Button disabled={parsedScheme.colors.length === 0} onClick={onApply} variant="contained">
                        Apply palette
                    </Button>
                </Stack>
            </DialogActions>
        </Dialog>
    );
}

export default CustomPaletteDialog;
