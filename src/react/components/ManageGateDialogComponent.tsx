import {
    Accordion,
    AccordionDetails,
    AccordionSummary,
    Divider,
    ListItemIcon,
    ListItemText,
    Menu,
    MenuItem,
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableRow,
    Tooltip,
    Typography,
} from "@mui/material";
import ExpandMoreIcon from "@mui/icons-material/ExpandMore";
import IconWithTooltip from "./IconWithTooltip";
import EditOutlinedIcon from "@mui/icons-material/EditOutlined";
import DeleteOutlineOutlinedIcon from "@mui/icons-material/DeleteOutlineOutlined";
import FileUploadOutlinedIcon from "@mui/icons-material/FileUploadOutlined";
import MoreVertIcon from "@mui/icons-material/MoreVert";
import { useCallback, useMemo, useState } from "react";
import { DEFAULT_GATE_COLOR, getRelevantGates } from "../gates/gateUtils";
import GateNameDialog from "./GateNameDialog";
import ConfirmDialog from "@/charts/dialogs/ConfirmDialog";
import { PenBoxIcon } from "lucide-react";
import { truncateGateLabel } from "../hooks/useGateLayers";
import { hexToRgb, rgbToHex } from "@/utilities/Utilities";
import type { Gate } from "../gates/types";
import { observer } from "mobx-react-lite";
import { useGateManager } from "../gates/useGateManager";

export type ManageGateDialogType = {
    open: boolean;
    onClose: () => void;
    xField?: string;
    yField?: string;
    onDelete: (gateId: string) => Promise<void>;
    onEdit: (gateId: string) => void;
    onRenameGate: (gateId: string, newName: string) => Promise<void>;
    onExportClick: (gateId: string) => void;
    onColorChange: (gateId: string, color: [number, number, number]) => Promise<void>;
    activeRegion?: string;
};

const ManageGateDialogContent = observer(({
    onClose,
    xField,
    yField,
    onDelete,
    onEdit,
    onRenameGate,
    onExportClick,
    onColorChange,
    activeRegion,
}: ManageGateDialogType) => {
    const gateManager = useGateManager();
    const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);
    const [selectedId, setSelectedId] = useState<string | null>(null);
    const [renameGateId, setRenameGateId] = useState<string | null>(null);
    const [deleteGateId, setDeleteGateId] = useState<string | null>(null);

    const relevantGates = useMemo(
        () =>
            getRelevantGates({
                gates: gateManager.gatesArray,
                xField,
                yField,
                region: activeRegion,
            }),
        [gateManager.gatesArray, xField, yField, activeRegion],
    );

    const regionSpecificGates = useMemo(() => 
        relevantGates.filter((gate) => gate.region === activeRegion), 
    [relevantGates, activeRegion]);

    const globalGates = useMemo(() => 
        relevantGates.filter((gate) => gate.region === undefined), 
    [relevantGates]);

    const handleMenuOpen = useCallback((event: React.MouseEvent<HTMLButtonElement>, gateId: string) => {
        event.stopPropagation();
        setSelectedId(gateId);
        setAnchorEl(event.currentTarget);
    }, []);

    const handleMenuClose = useCallback(() => {
        setAnchorEl(null);
        setSelectedId(null);
    }, []);

    const onRenameClick = useCallback((gateId: string) => {
        if (!gateId) return;
        setRenameGateId(gateId);
    }, []);

    const onDeleteClick = useCallback((gateId: string) => {
        if (!gateId) return;
        setDeleteGateId(gateId);
    }, []);

    const renderGateTable = useCallback(
        (gates: Gate[]) => (
            <Table
                sx={{
                    mb: 1,
                    width: "100%",
                    "& .MuiTableCell-head": {
                        fontWeight: 600,
                    },
                }}
            >
                <TableHead>
                    <TableRow>
                        <TableCell>
                            <Typography variant="button" color="text.secondary">
                                Name
                            </Typography>
                        </TableCell>
                        <TableCell sx={{ width: 48, minWidth: 48 }}>
                            <Typography variant="button" color="text.secondary">
                                Color
                            </Typography>
                        </TableCell>
                        <TableCell align="right" sx={{ width: 48, minWidth: 48 }}>
                            <Typography variant="button" color="text.secondary">
                                Actions
                            </Typography>
                        </TableCell>
                    </TableRow>
                </TableHead>
                <TableBody>
                    {gates.length === 0 ? (
                        <TableRow>
                            <TableCell colSpan={3} align="center" sx={{ py: 2 }}>
                                <Typography color="text.secondary">No gates in this section.</Typography>
                            </TableCell>
                        </TableRow>
                    ) : (
                        gates.map((gate) => (
                            <TableRow key={gate.id}>
                                <TableCell>
                                    <Tooltip title={gate.name}>
                                        <Typography variant="subtitle1">
                                            {truncateGateLabel(gate.name, 30)}
                                        </Typography>
                                    </Tooltip>
                                </TableCell>
                                <TableCell sx={{ width: 48, minWidth: 48, textAlign: "center" }}>
                                    <input
                                        type="color"
                                        value={rgbToHex(gate.color ?? DEFAULT_GATE_COLOR)}
                                        onChange={async (e) =>
                                            await onColorChange(gate.id, hexToRgb(e.target.value))
                                        }
                                        style={{
                                            width: 25,
                                            height: 25,
                                            padding: 0,
                                            borderRadius: 4,
                                            cursor: "pointer",
                                        }}
                                        aria-label={`Color for ${gate.name}`}
                                    />
                                </TableCell>
                                <TableCell align="right" sx={{ width: 48, minWidth: 48 }}>
                                    <IconWithTooltip
                                        tooltipText="More Actions"
                                        onClick={(e) => handleMenuOpen(e, gate.id)}
                                    >
                                        <MoreVertIcon />
                                    </IconWithTooltip>
                                </TableCell>
                            </TableRow>
                        ))
                    )}
                </TableBody>
            </Table>
        ),
        [onColorChange, handleMenuOpen],
    );

    return (
        <>
            {relevantGates.length === 0 ? (
                <Typography color="textSecondary" sx={{ py: 2 }}>
                    No gates yet. Draw a selection on the chart and save it as a gate.
                </Typography>
            ) : (
                <>
                    {activeRegion && (
                        <Accordion disableGutters sx={{ display: "block" }} defaultExpanded>
                            <AccordionSummary
                                expandIcon={<ExpandMoreIcon />}
                                sx={{ minHeight: 56 }}
                            >
                                <Typography variant="h6">
                                    "{activeRegion}" gates
                                </Typography>
                            </AccordionSummary>
                            <AccordionDetails sx={{ display: "block", pt: 0 }}>
                                {renderGateTable(regionSpecificGates)}
                            </AccordionDetails>
                        </Accordion>
                    )}

                    <Accordion disableGutters sx={{ display: "block", mt: 2 }} defaultExpanded>
                        <AccordionSummary
                            expandIcon={<ExpandMoreIcon />}
                            sx={{ minHeight: 56 }}
                        >
                            <Typography variant="h6">
                                Global gates
                            </Typography>
                        </AccordionSummary>
                        <AccordionDetails sx={{ display: "block", pt: 0 }}>
                            {renderGateTable(globalGates)}
                        </AccordionDetails>
                    </Accordion>
                </>
            )}
            {selectedId && (
                <Menu
                    anchorEl={anchorEl}
                    open={Boolean(anchorEl)}
                    onClose={handleMenuClose}
                    onClick={(e) => e.stopPropagation()}
                >
                    <MenuItem
                        onClick={() => {
                            onEdit(selectedId);
                            handleMenuClose();
                            onClose();
                        }}
                    >
                        <ListItemIcon>
                            <EditOutlinedIcon />
                        </ListItemIcon>
                        <ListItemText>Edit Geometry</ListItemText>
                    </MenuItem>
                    <MenuItem
                        onClick={() => {
                            onRenameClick(selectedId);
                            handleMenuClose();
                        }}
                    >
                        <ListItemIcon>
                            <PenBoxIcon size={20} />
                        </ListItemIcon>
                        <ListItemText>Rename</ListItemText>
                    </MenuItem>
                    <Divider />
                    <MenuItem
                        onClick={() => {
                            onExportClick(selectedId);
                            handleMenuClose();
                        }}
                    >
                        <ListItemIcon>
                            <FileUploadOutlinedIcon />
                        </ListItemIcon>
                        <ListItemText>Export (as GeoJSON)</ListItemText>
                    </MenuItem>
                    <Divider />
                    <MenuItem
                        onClick={() => {
                            onDeleteClick(selectedId);
                            handleMenuClose();
                        }}
                        sx={{ color: "var(--icon_color_error)" }}
                    >
                        <ListItemIcon sx={{ color: "var(--icon_color_error)" }}>
                            <DeleteOutlineOutlinedIcon />
                        </ListItemIcon>
                        <ListItemText>Delete</ListItemText>
                    </MenuItem>
                </Menu>
            )}

            {renameGateId && (
                <GateNameDialog
                    open={true}
                    onClose={() => setRenameGateId(null)}
                    onSaveGate={async (gateName) => {
                        await onRenameGate(renameGateId, gateName);
                    }}
                    name={relevantGates.find((g) => g.id === renameGateId)?.name ?? ""}
                    isEditName
                />
            )}

            {deleteGateId && (
                <ConfirmDialog
                    open={true}
                    title="Delete Gate"
                    message="Are you sure you want to delete the gate?"
                    onClose={() => setDeleteGateId(null)}
                    onConfirm={async () => {
                        await onDelete(deleteGateId);
                        setDeleteGateId(null);
                    }}
                    isDelete
                />
            )}
        </>
    );
});

export default ManageGateDialogContent;
