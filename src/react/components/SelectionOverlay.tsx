import { ButtonGroup, IconButton } from "@mui/material"; //check tree-shaking...
import PanToolOutlinedIcon from '@mui/icons-material/PanToolOutlined';
import PhotoSizeSelectSmallOutlinedIcon from '@mui/icons-material/PhotoSizeSelectSmallOutlined';
import PolylineOutlinedIcon from '@mui/icons-material/PolylineOutlined';
import EditOutlinedIcon from '@mui/icons-material/EditOutlined';
import { useMemo, useState } from "react";
import zIndex from "@mui/material/styles/zIndex";

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
    'lasso': {
        name: 'Lasso',
        ToolIcon: EditOutlinedIcon
    },
} as const;

type Tool = typeof Tools[keyof typeof Tools]['name'];
const ToolArray = Object.values(Tools);

function RectangleEditor({toolActive = false}) {
    // how shall we represent coordinates?
    const [start, setStart] = useState<[number, number]>([-100,-100]);
    const [end, setEnd] = useState<[number, number]>([0,0]);
    const [dragging, setDragging] = useState(false);

    return (
    <>
    <div style={{
        position: 'fixed',
        left: start[0],
        top: start[1],
        width: end[0] - start[0],
        height: end[1] - start[1],
        border: '1px solid white', //todo: make this a theme colour
        backgroundColor: 'rgba(255,255,255,0.1)',
        pointerEvents: 'none',
    }}/>
    <div style={{
        position: 'absolute',
        width: '100%',
        height: '100%',
    }} 
    onMouseDown={(e) => {
        if (!toolActive) return;
        setStart([e.clientX, e.clientY]);
        setEnd([e.clientX, e.clientY]);
        setDragging(true);
        window.addEventListener('mouseup', () => {
            setDragging(false);
        });
    }
    } onMouseMove={(e) => {
        if (dragging) {
            setEnd([e.clientX, e.clientY]);
        }
    }} onMouseUp={(e) => {
        setEnd([e.clientX, e.clientY]);
        setDragging(false);
    }} />
    </>);
}


export default function SelectionOverlay() {
    const [selectedTool, setSelectedTool] = useState<Tool>('Pan');
    // add a row of buttons to the top of the chart
    // rectangle, circle, polygon, lasso, magic wand, etc. 
    // (thanks copilot, that may be over-ambitious)
    // It would be good to have a poly-line tool with draggable points, though.
    // Also a brush tool for painting on masks with variable radius.
    const toolButtons = useMemo(() => {
        return ToolArray.map(({name, ToolIcon}) => {
            return <IconButton key={name} onClick={() => setSelectedTool(name)}><ToolIcon /></IconButton>
        });
    }, []);

    // state: { selectedTool: 'rectangle' | 'circle' | 'polygon' | 'lasso' | 'magic wand' | 'none' }
    // interaction phases... maybe revert back to pan after making selection
    // - but there should be interaction with drag handles...
    // add or remove from selection...
    // -> transform into Deck coordinates...
    // later: selection layers...
    return (
        <>
        <ButtonGroup variant="contained" aria-label="choose tool for manipulating view or selection" style={{zIndex: 2}}>
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
            onMouseDown={(e) => {
                
            }}
            >
                <RectangleEditor toolActive={selectedTool==='Rectangle'} />
        </div>
        </>
        )
}