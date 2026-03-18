import { ButtonGroup, Divider } from "@mui/material"; //check tree-shaking...
import PanToolOutlinedIcon from "@mui/icons-material/PanToolOutlined";
import PolylineOutlinedIcon from "@mui/icons-material/PolylineOutlined";
import PhotoSizeSelectSmallIcon from '@mui/icons-material/PhotoSizeSelectSmall';
import DoneOutlinedIcon from '@mui/icons-material/DoneOutlined';
import CloseOutlinedIcon from '@mui/icons-material/CloseOutlined';
import AddCircleOutlinedIcon from '@mui/icons-material/AddCircleOutlined';
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
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
import useGateActions from "../hooks/useGateActions";
import { useChartID, useParamColumns, useTheme } from "../hooks";
import IconWithTooltip from "./IconWithTooltip";
import { getEmptyFeatureCollection } from "../deck_state";
import { TuneOutlined } from "@mui/icons-material";
import { LassoIcon, SplineIcon } from "lucide-react";
import ManageGateDialogWrapper from "./ManageGateDialogWrapper";

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
        ToolIcon: PhotoSizeSelectSmallIcon,
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
    const chartId = useChartID();
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
    const manageGateByChartRef = useRef<Map<string, ManageGateDialogWrapper>>(new Map());
    const theme = useTheme();

    const {
        onDeleteGate,
        onExportClick,
        onRenameGate,
        onSaveGate,
        onEditGate,
        onConfirmEditGate,
        onCancelEditGate,
        onColorChange,
    } = useGateActions();

    useEffect(() => {
        // close the open dialog when unmounting the component
        return () => {
          const dlg = manageGateByChartRef.current.get(chartId);
          if (dlg) {
            setTimeout(() => dlg.close(), 0);
          }
          manageGateByChartRef.current.delete(chartId);
        };
      }, [chartId]);

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

    const relevantGates = useMemo(() => {
        if (!cx || !cy) return [];
        return gateManager.gatesArray.filter(
            (gate) => gate.columns[0] === cx.field && gate.columns[1] === cy.field
        )
    }
    , [gateManager.gatesArray, cx, cy]);

    const onManageGatesClick = useCallback(() => {
        console.log("chart id: ", chartId, manageGateByChartRef.current);
        if (manageGateByChartRef.current.get(chartId)) {
            // a dialog already exists
            return;
        }
        
        // Create a new manage gates dialog for a chart
        const dialog = new ManageGateDialogWrapper({
            gatesArray: relevantGates,
            onDelete: onDeleteGate,
            onColorChange,
            onRenameGate,
            onExportClick,
            onEdit: onEditGate,
            onClose: () => {
                // remove the entry of the dialog
                manageGateByChartRef.current.delete(chartId);
            }
        });

        // create a new entry for the dialog for a chart
        manageGateByChartRef.current.set(chartId, dialog);
    }, [
        relevantGates, 
        onDeleteGate, 
        onEditGate, 
        onColorChange, 
        onRenameGate, 
        onExportClick, 
        chartId
    ]);

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
                    onClick={onManageGatesClick}
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
                        tooltipText={"Save selection as Gate"}
                        onClick={() => setGateDialogOpen(true)}
                        iconButtonProps={{
                            "aria-label": "Save Gate",
                        }}
                    >
                        <AddCircleOutlinedIcon />
                    </IconWithTooltip>
                )}
                {editingGateId && (
                    <>
                        <IconWithTooltip
                            tooltipText={"Confirm"}
                            onClick={() => 
                                { onConfirmEditGate().catch(console.error) }
                            }
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
        </>
    );
});
