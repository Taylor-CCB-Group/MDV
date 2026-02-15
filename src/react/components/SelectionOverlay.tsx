import { ButtonGroup, Divider, IconButton, Tooltip } from "@mui/material"; //check tree-shaking...
import PanToolOutlinedIcon from "@mui/icons-material/PanToolOutlined";
import RectangleOutlinedIcon from "@mui/icons-material/RectangleOutlined";
import PolylineOutlinedIcon from "@mui/icons-material/PolylineOutlined";
import EditOutlinedIcon from "@mui/icons-material/EditOutlined";
import ControlCameraOutlinedIcon from "@mui/icons-material/ControlCameraOutlined";
import LayersIcon from '@mui/icons-material/Layers';
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
import { useTheme } from "../hooks";


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
const Tools = {
    pan: {
        name: "Pan",
        ToolIcon: PanToolOutlinedIcon,
        mode: TranslateModeEx
    },
    rectangle: {
        name: "Rectangle",
        ToolIcon: RectangleOutlinedIcon,
        mode: RectangleMode
    },
    // todo: add these back in once we have deck EditableGeoJsonLayer etc in place...
    polygon: {
        name: "Polygon",
        ToolIcon: PolylineOutlinedIcon,
        mode: PolygonMode
    },
    freehand: {
        name: "Freehand",
        ToolIcon: EditOutlinedIcon,
        mode: DrawPolygonByDraggingMode
    },
    transform: {
        name: "Modify",
        ToolIcon: ControlCameraOutlinedIcon,
        mode: EditMode
    },
    // measure: {
    //     name: "Measure",
    //     ToolIcon: StraightenIcon,
    //     mode: DrawLineStringMode
    // },
} as const;
type ToolIcon = typeof Tools[keyof typeof Tools]["ToolIcon"];
type Tool = (typeof Tools)[keyof typeof Tools]["name"];
const ToolArray = Object.values(Tools);

type ToolButtonProps = {
    name: Tool;
    ToolIcon: ToolIcon;
    selectedTool: Tool;
    setSelectedTool: (tool: Tool) => void;
};
// todo: fix the colors based on the themes
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
        <Tooltip title={name}>
            <IconButton
                sx={style}
                onClick={() => setSelectedTool(name)}
                aria-label={name}
            >
                <ToolIcon />
            </IconButton>
        </Tooltip>
    );
});

export default observer(function SelectionOverlay() {
    const { selectionProps } = useSpatialLayers();
    const { setSelectionMode, selectionFeatureCollection } = selectionProps;
    const gateManager = useGateManager();
    const [selectedTool, setSelectedToolX] = useState<Tool>("Pan");
    const [gateDialogOpen, setGateDialogOpen] = useState(false);
    const [manageGateDialogOpen, setManageGateDialogOpen] = useState(false);
    const theme = useTheme();

    const {
        onDeleteGate,
        onExportClick,
        onRenameGate,
        onSaveGate,
    } = useGateActions();

    const setSelectedTool = useCallback((tool: Tool) => {
        // pending refactor
        const mode = Object.values(Tools).find((t) => t.name === tool)?.mode;
        if (!mode) {
            console.error("no mode found for tool", tool);
            return;
        }
        //same composite mode order doesn't work for all tools, so making `mode()` be more explicit for each
        //setSelectionMode(new CompositeMode([new mode(), new TranslateModeEx()]));
        setSelectionMode(new mode());
        setSelectedToolX(tool);
    }, [setSelectionMode]);
    // add a row of buttons to the top of the chart
    // rectangle, circle, polygon, lasso, magic wand, etc.
    // (thanks copilot, that may be over-ambitious)
    // It would be good to have a poly-line tool with draggable points, though.
    // Also a brush tool for painting on masks with variable radius.
    const toolButtons = useMemo(() => {
        return ToolArray.map(({ name, ToolIcon, mode }) => (
            <ToolButton
                key={name}
                name={name}
                ToolIcon={ToolIcon}
                selectedTool={selectedTool}
                setSelectedTool={setSelectedTool}
            />
        ));
    }, [selectedTool, setSelectedTool]);

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
                <Tooltip title="Manage gates">
                    <IconButton
                        onClick={() => setManageGateDialogOpen(true)}
                        aria-label="Manage Gates"
                        sx={{
                            color: "var(--text_color)",
                            backgroundColor: "transparent",
                            borderRadius: 10,
                            zIndex: 2,
                            ml: 1,
                        }}
                    >
                        <LayersIcon />
                    </IconButton>
                </Tooltip>
                {hasSelection && (
                    <Tooltip title="Save selection as gate">
                        <IconButton
                            onClick={() => setGateDialogOpen(true)}
                            aria-label="Save as Gate"
                            sx={{
                                color: "var(--text_color)"
                            }}
                        >
                            <AddOutlinedIcon />
                        </IconButton>
                    </Tooltip>
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
                gatesArray={gateManager.gatesArray}
                onDelete={onDeleteGate}
                onRenameGate={onRenameGate}
                onExportClick={onExportClick}
            />
        </>
    );
});
