import {
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
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
    Typography,
} from "@mui/material";
import IconWithTooltip from "./IconWithTooltip";
import EditOutlinedIcon from "@mui/icons-material/EditOutlined";
import DeleteOutlineOutlinedIcon from "@mui/icons-material/DeleteOutlineOutlined";
import FileUploadOutlinedIcon from "@mui/icons-material/FileUploadOutlined";
import MoreVertIcon from "@mui/icons-material/MoreVert";
import { useCallback, useState } from "react";
import type { Gate } from "../gates/types";
import GateNameDialog from "./GateNameDialog";
import ConfirmDialog from "@/charts/dialogs/ConfirmDialog";
import { PenBoxIcon } from "lucide-react";
import { DialogCloseIconButton } from "@/catalog/ProjectRenameModal";

export type ManageGateDialogType = {
    open: boolean;
    onClose: () => void;
    gatesArray: Gate[];
    onDelete: (gateId: string) => void;
    onEdit: (gateId: string) => void;
    onRenameGate: (gateId: string, newName: string) => void;
    onExportClick: (gateId: string) => void;
};

const ManageGateDialog = ({
    open,
    onClose,
    gatesArray,
    onDelete,
    onEdit,
    onRenameGate,
    onExportClick,
}: ManageGateDialogType) => {
    const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);
    const [selectedId, setSelectedId] = useState<string | null>(null);
    const [renameGateId, setRenameGateId] = useState<string | null>(null);
    const [deleteGateId, setDeleteGateId] = useState<string | null>(null);

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

    return (
        <>
            <Dialog open={open} onClose={onClose} fullWidth maxWidth="xs">
                <DialogTitle>
                    Manage Gates
                    <DialogCloseIconButton onClose={onClose} />
                </DialogTitle>
                <DialogContent dividers>
                    <Table
                        sx={{
                            mb: 3,
                            "& .MuiTableCell-head": {
                                fontWeight: 600,
                            },
                        }}
                    >
                        <TableBody>
                            {gatesArray.length === 0 ? (
                                <TableRow>
                                <TableCell colSpan={2} align="center" sx={{ py: 4 }}>
                                    <Typography color="textSecondary">
                                        No gates yet. Draw a selection on the chart and save it as a gate.
                                    </Typography>
                                </TableCell>
                            </TableRow>
                            ) : ( 
                                gatesArray.map((gate) => (
                                <TableRow key={gate.id}>
                                    <TableCell>
                                        <Typography variant="subtitle1">{gate.name}</Typography>
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
                                    onRenameClick(selectedId)
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
                </DialogContent>
            </Dialog>
            {renameGateId && (
                <GateNameDialog
                    open={true}
                    onClose={() => {
                        setRenameGateId(null);
                    }}
                    onSaveGate={(gateName) => {
                        onRenameGate(renameGateId, gateName);
                    }}
                    name={gatesArray.find((g) => g.id === renameGateId)?.name ?? ""}
                />
            )}

            {deleteGateId && (
                <ConfirmDialog
                    open={true}
                    title="Delete Gate"
                    message="Are you sure you want to delete the gate?"
                    onClose={() => {
                        setDeleteGateId(null);
                    }}
                    onConfirm={() => {
                        onDelete(deleteGateId);
                        setDeleteGateId(null);
                    }}
                    isDelete
                />
            )}
        </>
    );
};

export default ManageGateDialog;
