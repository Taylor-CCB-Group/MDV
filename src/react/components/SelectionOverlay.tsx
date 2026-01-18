import { ButtonGroup, IconButton, Tooltip } from "@mui/material"; //check tree-shaking...
import PanToolOutlinedIcon from "@mui/icons-material/PanToolOutlined";
import PhotoSizeSelectSmallOutlinedIcon from "@mui/icons-material/PhotoSizeSelectSmallOutlined";
import PolylineOutlinedIcon from "@mui/icons-material/PolylineOutlined";
import EditOutlinedIcon from "@mui/icons-material/EditOutlined";
import ControlCameraOutlinedIcon from "@mui/icons-material/ControlCameraOutlined";
import BookmarkAddOutlinedIcon from "@mui/icons-material/BookmarkAddOutlined";
import { useCallback, useMemo, useState } from "react";
import type { useScatterplotLayer } from "../scatter_state";
import { useSpatialLayers } from "../spatial_context";
import type RangeDimension from "../../datastore/RangeDimension";
import { observer } from "mobx-react-lite";
import GateNameDialog from "./GateNameDialog";
import { useGateStore } from "../gates/useGateStore";
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
import type { Gate } from "../gates/types";
import { generateGateId } from "../gates/gateUtils";
import { useParamColumns } from "../hooks";
import { action } from "mobx";
import { useChart } from "../context";
import { getEmptyFeatureCollection } from "../deck_state";
import type { DeckScatterConfig } from "./DeckScatterReactWrapper";

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
        ToolIcon: PhotoSizeSelectSmallOutlinedIcon,
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
const ToolButton = ({ name, ToolIcon, selectedTool, setSelectedTool }: ToolButtonProps) => {
    const style = useMemo(
        () => ({
            backgroundColor:
                selectedTool === name ? "rgba(255,255,255,0.2)" : "transparent",
            borderRadius: 10,
            zIndex: 2,
        }),
        [selectedTool, name],
    );
    return (
        <Tooltip title={name}>
            <IconButton
                style={style}
                onClick={() => setSelectedTool(name)}
                aria-label={name}
            >
                <ToolIcon />
            </IconButton>
        </Tooltip>
    );
};

export default observer(function SelectionOverlay() {
    const { selectionProps } = useSpatialLayers();
    const { setSelectionMode, selectionFeatureCollection } = selectionProps;
    const gateStore = useGateStore();
    const paramColumns = useParamColumns();
    const chart = useChart<DeckScatterConfig>();
    const [selectedTool, setSelectedToolX] = useState<Tool>("Pan");
    const [gateDialogOpen, setGateDialogOpen] = useState(false);
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

    const onSaveGate = useCallback((gateName: string) => {
        if (gateStore.hasGateName(gateName)) {
            throw new Error('A gate with this name already exists');
        }
        
        // Get X and Y columns
        const [xCol, yCol] = paramColumns.slice(0, 2);
        if (!xCol || !yCol) {
            throw new Error('Chart must have X and Y axes');
        }
        
        // Create gate
        const gate: Gate = {
            id: generateGateId(),
            name: gateName.trim(),
            geometry: selectionFeatureCollection,
            columns: [xCol.field, yCol.field],
            createdAt: Date.now()
        };
        
        // Add to gate store
        gateStore.addGate(gate);
        
        // Clear selection
        action(() => {
            chart.config.selectionFeatureCollection = getEmptyFeatureCollection();
        })();
    }, [gateStore, paramColumns, selectionFeatureCollection, chart.config]);
    
    const hasSelection = useMemo(() => selectionFeatureCollection.features.length > 0, [selectionFeatureCollection.features.length]);
    
    return (
        <>
            <ButtonGroup
                variant="contained"
                aria-label="choose tool for manipulating view or selection"
                //moving this to the top right corner and absolute to avoid interfering with axes
                className="z-[2] p-2 absolute top-0 right-0"
                // style={{zIndex: 2, padding: '0.3em'}}
            >
                {toolButtons}
                {hasSelection && (
                    <Tooltip title="Save selection as gate">
                        <IconButton
                            onClick={() => setGateDialogOpen(true)}
                            aria-label="Save as Gate"
                            style={{ backgroundColor: 'rgba(76, 175, 80, 0.2)' }}
                        >
                            <BookmarkAddOutlinedIcon />
                        </IconButton>
                    </Tooltip>
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
