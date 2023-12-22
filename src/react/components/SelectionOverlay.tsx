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

export default function SelectionOverlay() {
    const [selectedTool, setSelectedTool] = useState<Tool>('Pan');
    // add a row of buttons to the top of the chart
    // rectangle, circle, polygon, lasso, magic wand, etc. 
    // (thanks copilot, that may be over-ambitious)
    // It would be good to have a poly-line tool with draggable points, though.
    const toolButtons = useMemo(() => {
        return ToolArray.map(({name, ToolIcon}) => {
            return <IconButton key={name} onClick={() => setSelectedTool(name)}><ToolIcon /></IconButton>
        });
    }, []);

    // state: { selectedTool: 'rectangle' | 'circle' | 'polygon' | 'lasso' | 'magic wand' | 'none' }
    // interaction phases...
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
            // backgroundColor: 'red',
            // opacity: 0.5,
            width: '100%',
            height: '100%',
            zIndex: 1,
            pointerEvents: selectedTool === 'Pan' ? 'none' : 'auto'
            }}>
        </div>
        </>
        )
}