import {
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
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
import BorderColorOutlinedIcon from "@mui/icons-material/BorderColorOutlined";
import DeleteOutlineOutlinedIcon from "@mui/icons-material/DeleteOutlineOutlined";
import FileUploadOutlinedIcon from "@mui/icons-material/FileUploadOutlined";
import MoreVertIcon from "@mui/icons-material/MoreVert";
import { useCallback, useState } from "react";
import type { Gate } from "../gates/types";
import GateNameDialog from "./GateNameDialog";
import ConfirmDialog from "@/charts/dialogs/ConfirmDialog";

export type ManageGateDialogType = {
    open: boolean;
    onClose: () => void;
    gatesArray: Gate[];
    onDelete: (gateId: string) => void;
    onRenameGate: (gateId: string, newName: string) => void;
    onExportClick: (gateId: string) => void;
};

const ManageGateDialog = ({
    open,
    onClose,
    gatesArray,
    onDelete,
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
            <Dialog open={open} onClose={onClose} fullWidth maxWidth="sm">
                <DialogContent dividers>
                    <Table
                        sx={{
                            mb: 3,
                            "& .MuiTableCell-head": {
                                fontWeight: 600,
                            },
                        }}
                    >
                        <TableHead>
                            <TableRow>
                                <TableCell>
                                    <Typography variant="button" color="textSecondary">
                                        Name
                                    </Typography>
                                </TableCell>
                                <TableCell align="right" sx={{ width: 48, minWidth: 48 }}>
                                    <Typography variant="button" color="textSecondary">
                                        Rename
                                    </Typography>
                                </TableCell>
                                <TableCell align="right" sx={{ width: 48, minWidth: 48 }}>
                                    <Typography variant="button" color="textSecondary">
                                        Delete
                                    </Typography>
                                </TableCell>
                                <TableCell align="right" sx={{ width: 48, minWidth: 48 }}>
                                    <Typography variant="button" color="textSecondary">
                                        Actions
                                    </Typography>
                                </TableCell>
                            </TableRow>
                        </TableHead>
                        <TableBody>
                            {gatesArray.map((gate) => (
                                <TableRow key={gate.id}>
                                    <TableCell>
                                        <Typography variant="subtitle1">{gate.name}</Typography>
                                    </TableCell>
                                    <TableCell align="right" sx={{ width: 48, minWidth: 48 }}>
                                        {/* todo: add this onClick*/}
                                        <IconWithTooltip
                                            tooltipText="Rename Gate"
                                            onClick={() => onRenameClick(gate.id)}
                                        >
                                            <BorderColorOutlinedIcon />
                                        </IconWithTooltip>
                                    </TableCell>
                                    <TableCell align="right" sx={{ width: 48, minWidth: 48 }}>
                                        <IconWithTooltip
                                            tooltipText="Delete Gate"
                                            onClick={() => onDeleteClick(gate.id)}
                                        >
                                            <DeleteOutlineOutlinedIcon />
                                        </IconWithTooltip>
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
                            ))}
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
                                    onExportClick(selectedId);
                                    handleMenuClose();
                                }}
                            >
                                <ListItemIcon>
                                    <FileUploadOutlinedIcon />
                                </ListItemIcon>
                                <ListItemText>Export (as GeoJSON)</ListItemText>
                            </MenuItem>
                        </Menu>
                    )}
                </DialogContent>
                <DialogActions>
                    <Button onClick={onClose} color="error">
                        Close
                    </Button>
                </DialogActions>
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
                    }}
                />
            )}
        </>
    );
};

export default ManageGateDialog;
