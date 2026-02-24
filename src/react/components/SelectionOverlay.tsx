import { ButtonGroup, Divider } from "@mui/material"; //check tree-shaking...
import PanToolOutlinedIcon from "@mui/icons-material/PanToolOutlined";
import RectangleOutlinedIcon from "@mui/icons-material/RectangleOutlined";
import PolylineOutlinedIcon from "@mui/icons-material/PolylineOutlined";
import EditOutlinedIcon from "@mui/icons-material/EditOutlined";
import ControlCameraOutlinedIcon from "@mui/icons-material/ControlCameraOutlined";
import DoneOutlinedIcon from '@mui/icons-material/DoneOutlined';
import CloseOutlinedIcon from '@mui/icons-material/CloseOutlined';
import AddOutlinedIcon from '@mui/icons-material/AddOutlined';
import { useCallback, useMemo, useState } from "react";
import { useSpatialLayers } from "../spatial_context";
import { observer } from "mobx-react-lite";
import GateNameDialog from "./GateNameDialog";
import { useGateManager } from "../gates/useGateManager";
import {
    DrawPolygonMode,
    DrawPolygonByDraggingMode,
    ModifyMode,
    // DrawRectangleMode,
    // TransformMode,
    // TranslateMode,
    CompositeMode
} from '@deck.gl-community/editable-layers';
import TranslateModeEx from '../../editable-layers/deck-community-ish/translate-mode-exp';
import { DrawRectangleByDraggingMode } from "@/editable-layers/deck-community-ish/draw-rectangle-by-dragging-mode";
import ManageGateDialog from "./ManageGateDialog";
import useGateActions from "../hooks/useGateActions";
import { useParamColumns, useTheme } from "../hooks";
import IconWithTooltip from "./IconWithTooltip";
import { getEmptyFeatureCollection } from "../deck_state";
import { TuneOutlined } from "@mui/icons-material";
import { LassoIcon, PentagonIcon, SplineIcon, SplinePointerIcon, SquareIcon } from "lucide-react";

class EditMode extends CompositeMode {
    constructor() {
        super([
            new TranslateModeEx(), //works with cartesian / non-GIS coordinates
            new ModifyMode(),
        ]);
    }
}


class RectangleMode extends CompositeMode {
    constructor() {
        super([
            // we pass in modeProps to the layer, not the edit mode.
            // new DrawRectangleMode(), //with modeConfig: { dragToDraw: true }
            //our version, which we shouldn't need but looks/works a bit different as of now
            new DrawRectangleByDraggingMode(),
            new TranslateModeEx(),
        ]);
    }
}
class PolygonMode extends CompositeMode {
    constructor() {
        super([
            new DrawPolygonMode(),
            new TranslateModeEx(),
        ]);
    }
}
// todo - figure out weird conflict behaviour with this mode...
class FreehandMode extends CompositeMode {
    constructor() {
        super([
            new DrawPolygonByDraggingMode(),
            new TranslateModeEx(),
        ]);
    }
}
// material-ui icons, or font-awesome icons... or custom in some cases...
// mui icons are hefty, not sure about this...
export const Tools = {
    pan: {
        name: "Pan",
        ToolIcon: PanToolOutlinedIcon,
        mode: TranslateModeEx
    },
    rectangle: {
        name: "Rectangle",
        ToolIcon: SquareIcon,
        mode: RectangleMode
    },
    // todo: add these back in once we have deck EditableGeoJsonLayer etc in place...
    polygon: {
        name: "Polygon",
        ToolIcon: PentagonIcon,
        mode: PolygonMode
    },
    freehand: {
        name: "Freehand",
        ToolIcon: LassoIcon,
        mode: DrawPolygonByDraggingMode
    },
    transform: {
        name: "Modify",
        ToolIcon: SplineIcon,
        mode: EditMode
    },
    // measure: {
    //     name: "Measure",
    //     ToolIcon: StraightenIcon,
    //     mode: DrawLineStringMode
    // },
} as const;
type ToolIcon = typeof Tools[keyof typeof Tools]["ToolIcon"];
export type Tool = (typeof Tools)[keyof typeof Tools]["name"];
const ToolArray = Object.values(Tools);
const EDIT_MODE_TOOL_NAMES: Tool[] = ["Pan", "Modify"];

type ToolButtonProps = {
    name: Tool;
    ToolIcon: ToolIcon;
    selectedTool: Tool;
    setSelectedTool: (tool: Tool) => void;
};

const ToolButton = observer(({ name, ToolIcon, selectedTool, setSelectedTool }: ToolButtonProps) => {
    const theme = useTheme();
    const selectedBackgroundColor = theme === "dark" ? "rgba(255,255,255,0.2)" : "rgba(48, 46, 46, 0.18)";
    const style = useMemo(
        () => ({
            backgroundColor:
                selectedTool === name ? selectedBackgroundColor : "transparent",
            color: "var(--text_color)",
            borderRadius: 10,
            zIndex: 2,
        }),
        [selectedTool, name, selectedBackgroundColor],
    );
    return (
        <IconWithTooltip
            tooltipText={name}
            onClick={() => setSelectedTool(name)}
            iconButtonProps={{
                sx: style,
                "aria-label": name,
            }}
        >
            <ToolIcon size={20} />
        </IconWithTooltip>
    );
});

export default observer(function SelectionOverlay() {
    const { selectionProps } = useSpatialLayers();
    const { 
        setSelectionMode, 
        selectionFeatureCollection,
        editingGateId, 
        setSelectionFeatureCollection,
        selectedTool,
        setSelectedTool: setSelectedToolX
    } = selectionProps;
    const gateManager = useGateManager();
    const [cx, cy] = useParamColumns();
    
    const [gateDialogOpen, setGateDialogOpen] = useState(false);
    const [manageGateDialogOpen, setManageGateDialogOpen] = useState(false);
    const theme = useTheme();

    const {
        onDeleteGate,
        onExportClick,
        onRenameGate,
        onSaveGate,
        onEditGate,
        onConfirmEditGate,
        onCancelEditGate,
    } = useGateActions();

    const toolsToShow = useMemo(
        () =>
            editingGateId
                ? ToolArray.filter((t) => EDIT_MODE_TOOL_NAMES.includes(t.name))
                : ToolArray,
        [editingGateId],
    );

    const setSelectedTool = useCallback((tool: Tool) => {
        // pending refactor
        const mode = Object.values(Tools).find((t) => t.name === tool)?.mode;
        if (!mode) {
            console.error("no mode found for tool", tool);
            return;
        }

        if (!editingGateId) {
            setSelectionFeatureCollection(getEmptyFeatureCollection());
        }
        //same composite mode order doesn't work for all tools, so making `mode()` be more explicit for each
        //setSelectionMode(new CompositeMode([new mode(), new TranslateModeEx()]));
        setSelectionMode(new mode());
        setSelectedToolX(tool);
    }, [setSelectionMode, editingGateId, setSelectionFeatureCollection, setSelectedToolX]);
    // add a row of buttons to the top of the chart
    // rectangle, circle, polygon, lasso, magic wand, etc.
    // (thanks copilot, that may be over-ambitious)
    // It would be good to have a poly-line tool with draggable points, though.
    // Also a brush tool for painting on masks with variable radius.
    const toolButtons = useMemo(() => {
        return toolsToShow.map(({ name, ToolIcon }) => (
            <ToolButton
                key={name}
                name={name}
                ToolIcon={ToolIcon}
                selectedTool={selectedTool}
                setSelectedTool={setSelectedTool}
            />
        ));
    }, [selectedTool, setSelectedTool, toolsToShow]);

    const relevantGates = useMemo(() => 
        gateManager.gatesArray.filter(
            (gate) => gate.columns[0] === cx.field && gate.columns[1] === cy.field
        )
    , [gateManager.gatesArray, cx, cy]);

    const hasSelection = useMemo(() => selectionFeatureCollection.features.length > 0, [selectionFeatureCollection.features.length]);
    
    return (
        <>
            <ButtonGroup
                variant="contained"
                aria-label="choose tool for manipulating view or selection"
                //moving this to the top right corner and absolute to avoid interfering with axes
                className="z-[2] p-2 absolute top-0 right-0 opacity-90"
                sx={{
                    backgroundColor: theme === "dark" ? "rgba(37,37,37,0.92)" : "rgba(245,245,245,0.92)",
                    borderRadius: 0,
                    borderBottomLeftRadius: "5px",
                }}
            >
                {toolButtons}
                <Divider 
                    orientation="vertical" 
                    sx={{
                        color: "inherit",
                        width: "5px",
                        height: "35px",
                        padding: "2px",
                    }} 
                />
                <IconWithTooltip
                    tooltipText={"Manage Gates"}
                    onClick={() => setManageGateDialogOpen(true)}
                    iconButtonProps={{
                        sx: {
                            color: "var(--text_color)",
                            borderRadius: 10,
                            zIndex: 2,
                            ml: 1,
                        },
                        "aria-label": "Manage Gates",
                    }}
                >
                    <TuneOutlined />
                </IconWithTooltip>
                {hasSelection && !editingGateId && (
                    <IconWithTooltip
                        tooltipText={"Save selection as gate"}
                        onClick={() => setGateDialogOpen(true)}
                        iconButtonProps={{
                            "aria-label": "Save Gate",
                        }}
                    >
                        <AddOutlinedIcon color="success" />
                    </IconWithTooltip>
                )}
                {editingGateId && (
                    <>
                        <IconWithTooltip
                            tooltipText={"Confirm"}
                            onClick={onConfirmEditGate}
                            iconButtonProps={{
                                sx: {
                                    color: "green",
                                    borderRadius: 10,
                                    zIndex: 2,
                                },
                                "aria-label": "Confirm",
                            }}
                        >
                            <DoneOutlinedIcon />
                        </IconWithTooltip>
                        <IconWithTooltip
                            tooltipText={"Cancel"}
                            onClick={onCancelEditGate}
                            iconButtonProps={{
                                sx: {
                                    color: "red",
                                    borderRadius: 10,
                                    zIndex: 2,
                                },
                                "aria-label": "Cancel",
                            }}
                        >
                            <CloseOutlinedIcon />
                        </IconWithTooltip>
                    </>
                )}
            </ButtonGroup>
            
            <GateNameDialog
                open={gateDialogOpen}
                onClose={() => setGateDialogOpen(false)}
                onSaveGate={onSaveGate}
            />
            <ManageGateDialog
                open={manageGateDialogOpen}
                onClose={() => setManageGateDialogOpen(false)}
                // Passing gatesArray directly as there are no mutations to the array
                gatesArray={relevantGates}
                onDelete={onDeleteGate}
                onRenameGate={onRenameGate}
                onExportClick={onExportClick}
                onEdit={onEditGate}
            />
        </>
    );
});
