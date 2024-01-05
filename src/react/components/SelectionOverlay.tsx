import { ButtonGroup, IconButton, Tooltip } from "@mui/material"; //check tree-shaking...
import PanToolOutlinedIcon from '@mui/icons-material/PanToolOutlined';
import PhotoSizeSelectSmallOutlinedIcon from '@mui/icons-material/PhotoSizeSelectSmallOutlined';
import PolylineOutlinedIcon from '@mui/icons-material/PolylineOutlined';
import EditOutlinedIcon from '@mui/icons-material/EditOutlined';
import { useCallback, useEffect, useMemo, useState } from "react";
import { useViewerStore } from "./avivatorish/state";
import { ScatterplotLayer } from "deck.gl/typed";

// material-ui icons, or font-awesome icons... or custom in some cases...
// mui icons are hefty, not sure about this...
const Tools = {
    'pan': {
        name: 'Pan',
        ToolIcon: PanToolOutlinedIcon
    },
    'rectangle': {
        name: 'Rectangle',
        ToolIcon: PhotoSizeSelectSmallOutlinedIcon //todo something better...
    },
    'polygon': {
        name: 'Polygon',
        ToolIcon: PolylineOutlinedIcon
    },
    'freehand': {
        name: 'Freehand',
        ToolIcon: EditOutlinedIcon
    },
} as const;

type Tool = typeof Tools[keyof typeof Tools]['name'];
const ToolArray = Object.values(Tools);

function RectangleEditor({toolActive = false, scatterplotLayer} : {toolActive: boolean, scatterplotLayer: ScatterplotLayer}) {
    // how shall we represent coordinates?
    // - should be in Deck coordinates, we need to convert to/from screen coordinates
    // - we may still want to use screen coordinates locally, but there should be a store
    // with the Deck coordinates and also methods for actually doing the selection etc.
    // (may well be in a worker - or perhaps we can use the GPU for this)
    const [start, setStart] = useState<[number, number]>([0,0]);
    const [end, setEnd] = useState<[number, number]>([0,0]);
    const [dragging, setDragging] = useState(false);
    const handleMouseUp = useCallback(() => {
        setDragging(false);
        window.removeEventListener('mouseup', handleMouseUp);
    }, []);
    const unproject = useCallback((e: React.MouseEvent) => {
        const r = e.currentTarget.getBoundingClientRect();
        const x = e.clientX - r.left;
        const y = e.clientY - r.top;
        const p = scatterplotLayer.unproject([x, y]) as [number, number];
        return p;
    }, [scatterplotLayer]);
    const viewState = useViewerStore((state) => state.viewState); // for reactivity - still a frame behind...
    const [frameTrigger, setFrameTrigger] = useState(false);
    useEffect(() => {
        //non-trivial overhead here? must be a better way...
        const t = setTimeout(() => {
            setFrameTrigger(v => !v);
        }, 0);
        return () => clearTimeout(t);
    }, [viewState]);
    const min = [Math.min(start[0], end[0]), Math.min(start[1], end[1])];
    const max = [Math.max(start[0], end[0]), Math.max(start[1], end[1])];
    const screenStart = scatterplotLayer.project(min);
    const screenEnd = scatterplotLayer.project(max);

    return (
    <>
    <div style={{
        // position: dragging ? 'fixed' : 'relative',
        position: 'relative',
        left: screenStart[0],
        top: screenStart[1],
        width: screenEnd[0] - screenStart[0],
        height: screenEnd[1] - screenStart[1],
        border: '1px solid white', //todo: make this a theme colour
        backgroundColor: 'rgba(255,255,255,0.1)',
    }}
    />    
    <div style={{
        position: 'absolute',
        width: '100%',
        height: '100%',
        top: 0,
        left: 0,
    }} 
    onMouseDown={(e) => {
        if (!toolActive) return;
        const p = unproject(e);
        setStart(p);
        setEnd(p);
        setDragging(true);
        window.addEventListener('mouseup', handleMouseUp);
    }
    } onMouseMove={(e) => {
        if (!toolActive) return;
        if (dragging) {
            const p = unproject(e);
            setEnd(p);
        }
    }} onMouseUp={(e) => {
        if (!toolActive) return;
        const p = unproject(e);
        setEnd(p);
        setDragging(false);
    }} />
    </>);
}


export default function SelectionOverlay({scatterplotLayer} : {scatterplotLayer: ScatterplotLayer}) {
    //for now, passing down scatterplotLayer so we can use it to unproject screen coordinates, 
    //as we figure out how to knit this together...
    const [selectedTool, setSelectedTool] = useState<Tool>('Pan');
    // add a row of buttons to the top of the chart
    // rectangle, circle, polygon, lasso, magic wand, etc. 
    // (thanks copilot, that may be over-ambitious)
    // It would be good to have a poly-line tool with draggable points, though.
    // Also a brush tool for painting on masks with variable radius.
    const toolButtons = useMemo(() => {
        return ToolArray.map(({name, ToolIcon}) => {
            return <Tooltip title={name} key={name}><IconButton style={{
                backgroundColor: selectedTool === name ? 'rgba(255,255,255,0.2)' : 'transparent',
                borderRadius: 10,
                zIndex: 2,
            }} onClick={() => setSelectedTool(name)}
            aria-label={name}
            ><ToolIcon /></IconButton></Tooltip>
        });
    }, [selectedTool]);
    // state: { selectedTool: 'rectangle' | 'circle' | 'polygon' | 'lasso' | 'magic wand' | 'none' }
    // interaction phases... maybe revert back to pan after making selection
    // - but there should be interaction with drag handles...
    // add or remove from selection...
    // -> transform into Deck coordinates...
    // later: selection layers...
    const drawRect = scatterplotLayer && scatterplotLayer.internalState;
    return (
        <>
        <ButtonGroup variant="contained" aria-label="choose tool for manipulating view or selection" style={{zIndex: 2, padding: '0.3em'}}>
            {toolButtons}
        </ButtonGroup>
        <div style={{
            position: 'absolute', 
            top: 0,
            width: '100%',
            height: '100%',
            zIndex: 1,
            pointerEvents: selectedTool === 'Pan' ? 'none' : 'auto'
            }}
            onMouseUp={(e) => {
                setSelectedTool('Pan');
            }}
            onKeyDown={(e) => {
                //not working as of now... probably want to think more about keyboard shortcuts in general...
                if (e.key === 'Escape') {
                    setSelectedTool('Pan');
                }
                if (e.key === 'r') {
                    setSelectedTool('Rectangle');
                }
            }}
            >
                {drawRect && <RectangleEditor toolActive={selectedTool==='Rectangle'} scatterplotLayer={scatterplotLayer}/>}
        </div>
        </>
    )
}