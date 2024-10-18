import { ButtonGroup, IconButton, Tooltip } from "@mui/material"; //check tree-shaking...
import PanToolOutlinedIcon from "@mui/icons-material/PanToolOutlined";
import PhotoSizeSelectSmallOutlinedIcon from "@mui/icons-material/PhotoSizeSelectSmallOutlined";
import PolylineOutlinedIcon from "@mui/icons-material/PolylineOutlined";
import EditOutlinedIcon from "@mui/icons-material/EditOutlined";
import ControlCameraOutlinedIcon from "@mui/icons-material/ControlCameraOutlined";
import StraightenIcon from "@mui/icons-material/Straighten";
import { useCallback, useMemo, useRef, useState } from "react";
import { useMetadata } from "./avivatorish/state";
import { useRegionScale, type useScatterplotLayer } from "../scatter_state";
import { useMeasure, useSpatialLayers } from "../spatial_context";
import type RangeDimension from "../../datastore/RangeDimension";
import { observer } from "mobx-react-lite";
import { useChartDoc, useChartSize } from "../hooks";
import { sizeToMeters } from "./avivatorish/utils";
import clsx from "clsx";
import {
    DrawPolygonMode,
    DrawPolygonByDraggingMode,
    ModifyMode,
    DrawRectangleMode, //we need a different version of this that works by dragging
    DrawLineStringMode,
    // TransformMode,
    // TranslateMode,
    CompositeMode,
} from '@deck.gl-community/editable-layers';
import TranslateModeEx from '../../editable-layers/deck-community-ish/translate-mode-exp';
import { DrawRectangleByDraggingMode } from "@/editable-layers/deck-community-ish/draw-rectangle-by-dragging-mode";

class EditMode extends CompositeMode {
    constructor() {
        super([
            new TranslateModeEx(), //works with cartesian / non-GIS coordinates
            new ModifyMode(),
        ]);
    }
}

//might want to have EdtitMode for this, could be some confusion with this & previous transform mode?
class PanMode extends CompositeMode {
    constructor() {
        super([]);
    }
}

// material-ui icons, or font-awesome icons... or custom in some cases...
// mui icons are hefty, not sure about this...
const Tools = {
    pan: {
        name: "Pan",
        ToolIcon: PanToolOutlinedIcon,
        mode: PanMode
    },
    rectangle: {
        name: "Rectangle",
        ToolIcon: PhotoSizeSelectSmallOutlinedIcon, //todo something better...
        mode: DrawRectangleByDraggingMode
    },
    // todo: add these back in once we have deck EditableGeoJsonLayer etc in place...
    polygon: {
        name: "Polygon",
        ToolIcon: PolylineOutlinedIcon,
        mode: DrawPolygonMode
    },
    freehand: {
        name: "Freehand",
        ToolIcon: EditOutlinedIcon,
        mode: DrawPolygonByDraggingMode
    },
    transform: {
        name: "Transform",
        ToolIcon: ControlCameraOutlinedIcon,
        mode: EditMode
    },
    measure: {
        name: "Measure",
        ToolIcon: StraightenIcon,
        mode: DrawLineStringMode
    },
} as const;

type Tool = (typeof Tools)[keyof typeof Tools]["name"];
const ToolArray = Object.values(Tools);
type P = [number, number];
type EditorProps = {
    toolActive?: boolean;
    rangeDimension: RangeDimension;
} & ReturnType<typeof useScatterplotLayer>;

function MeasureTool({ scatterplotLayer, unproject, toolActive }: EditorProps) {
    // click to set start, click to set end, draw line between them.
    // UX for tools is generally in need of work...
    const scale = useRegionScale();
    const metadata = useMetadata();
    const physicalSize = metadata?.Pixels?.PhysicalSizeX || 1;
    // const physicalUnits = metadata?.Pixels?.PhysicalSizeXUnit || 'unknown units';
    const [w, h] = useChartSize();
    const [hasStarted, setHasStarted] = useState(false);
    const { startPixels, setStart, endPixels, setEnd } = useMeasure();
    const handleMouseMove = useCallback(
        (e: MouseEvent) => {
            setEnd(unproject(e));
        },
        [unproject, setEnd],
    );
    const doc = useChartDoc();
    const handleMouseUp = useCallback(
        (e: MouseEvent) => {
            handleMouseMove(e);
            doc.removeEventListener("mouseup", handleMouseUp);
            doc.removeEventListener("mousemove", handleMouseMove);
        },
        [handleMouseMove, doc],
    );
    if (!toolActive && !hasStarted) return null;
    // could we unproject into the image layer rather than scatterplotLayer?
    // const startPixels = unproject(start);
    // const endPixels = unproject(end);
    // todo make this a layer in deck.gl?
    console.warn("MeasureTool is in development and may not work as expected.");
    const start = scatterplotLayer.project(startPixels);
    const end = scatterplotLayer.project(endPixels);
    const length =
        (physicalSize *
            Math.sqrt(
                (endPixels[0] - startPixels[0]) ** 2 +
                    (endPixels[1] - startPixels[1]) ** 2,
            )) /
        scale;
    const f = Math.round;
    // <<< XXX this isn't working and I don't know why... values look right but css isn't applied >>>
    // const className = clsx("absolute", `top-[${f(end[1] + 10)}px]`, `left-[${f(end[0] + 10)}px]`); //possible HMR problem with rule-of-hooks?
    // console.log('measure tool className=', className);
    return (
        <>
            <div
                // todo make this a layer in deck.gl? or at least avoid style literal...
                // className={clsx("absolute", `top-[${f(end[1] + 10)}px]`, `left-[${f(end[0] + 10)}px]`)}
                // className={clsx("absolute", "top-[380px]", `left-[${f(end[0] + 10)}px]`)} //this works
                // className={className}

                style={{
                    position: "absolute",
                    width: "100%",
                    height: "100px",
                    display: hasStarted ? "block" : "none",
                    top: end[1] + 10,
                    left: end[0] + 10,
                }}
            >
                {sizeToMeters(length, "mm").toFixed(2)}mm
            </div>
            <svg
                viewBox={`0 0 ${w} ${h}`}
                className={clsx(
                    "absolute top-0 left-0 w-full h-full bg-transparent",
                    toolActive ? "pointer-events-auto" : "pointer-events-none",
                )}
                onMouseDown={(e) => {
                    if (!toolActive) return;
                    const p = unproject([e.clientX, e.clientY]);
                    setStart(p);
                    setEnd(p);
                    setHasStarted(true);
                    doc.addEventListener("mouseup", handleMouseUp);
                    doc.addEventListener("mousemove", handleMouseMove);
                }}
                onScroll={(e) => {
                    e.preventDefault();
                }}
            >
                {hasStarted && (
                    <>
                        <circle
                            cx={start[0]}
                            cy={start[1]}
                            r={3}
                            fill="white"
                        />
                        <circle cx={end[0]} cy={end[1]} r={3} fill="white" />
                        <line
                            x1={start[0]}
                            y1={start[1]}
                            x2={end[0]}
                            y2={end[1]}
                            stroke="white"
                            strokeWidth={1.5}
                        />
                    </>
                )}
            </svg>
        </>
    );
}

function TransformEditor({
    scatterplotLayer,
    modelMatrix,
    unproject,
}: EditorProps) {
    const pLastRef = useRef([Number.NaN, Number.NaN]);
    const handleMouseMove = useCallback(
        (e: MouseEvent) => {
            const p = unproject(e);
            const pLast = pLastRef.current;
            const dx = p[0] - pLast[0];
            const dy = p[1] - pLast[1];
            modelMatrix.translate([dx, dy, 0]);
            //this may redraw nice & fast without incurring react overhead, but now we're left with other problems i.e. rectangle editor not updating...
            scatterplotLayer.setNeedsRedraw();
        },
        [modelMatrix.translate, scatterplotLayer.setNeedsRedraw, unproject],
    );
    const doc = useChartDoc();
    const handleMouseUp = useCallback(
        (e: MouseEvent) => {
            doc.removeEventListener("mouseup", handleMouseUp);
            doc.removeEventListener("mousemove", handleMouseMove);
        },
        [doc, handleMouseMove],
    );

    return (
        <>
            <div
                className="absolute top-0 left-0 w-full h-full"
                onMouseDown={(e) => {
                    const p = unproject(e);
                    pLastRef.current = p;
                    doc.addEventListener("mouseup", handleMouseUp);
                    doc.addEventListener("mousemove", handleMouseMove);
                }}
            />
        </>
    );
}

const ToolButton = ({ name, ToolIcon, selectedTool, setSelectedTool }) => {
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
    const { setSelectionMode } = selectionProps;
    const [selectedTool, setSelectedToolX] = useState<Tool>("Pan");
    const setSelectedTool = useCallback((tool: Tool) => {
        // pending refactor
        const mode = Object.values(Tools).find((t) => t.name === tool)?.mode;
        if (!mode) {
            console.error("no mode found for tool", tool);
            return;
        }
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
    }, [selectedTool]);
    return (
        <>
            <ButtonGroup
                variant="contained"
                aria-label="choose tool for manipulating view or selection"
                className="z-[2] p-2"
                // style={{zIndex: 2, padding: '0.3em'}}
            >
                {toolButtons}
            </ButtonGroup>
        </>
    );
});
