import { ButtonGroup, IconButton, Tooltip } from "@mui/material"; //check tree-shaking...
import PanToolOutlinedIcon from '@mui/icons-material/PanToolOutlined';
import PhotoSizeSelectSmallOutlinedIcon from '@mui/icons-material/PhotoSizeSelectSmallOutlined';
import PolylineOutlinedIcon from '@mui/icons-material/PolylineOutlined';
import EditOutlinedIcon from '@mui/icons-material/EditOutlined';
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useViewerStore } from "./avivatorish/state";
import { ScatterplotLayer } from "deck.gl/typed";
import { useRegionScale } from "../scatter_state";
import { useChart } from "../context";
import RangeDimension from "../../datastore/RangeDimension";
import { useChartID } from "../hooks";
import { observer } from "mobx-react-lite";

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
type P = [number, number];
type EditorProps = { toolActive: boolean, scatterplotLayer: ScatterplotLayer, rangeDimension: RangeDimension };
function RectangleEditor({toolActive = false, scatterplotLayer, rangeDimension} : EditorProps) {
    // how shall we represent coordinates?
    // - should be in Deck coordinates, we need to convert to/from screen coordinates
    // - we may still want to use screen coordinates locally, but there should be a store
    // with the Deck coordinates and also methods for actually doing the selection etc.
    // (may well be in a worker - or perhaps we can use the GPU for this)
    // - should be a RangeDimension, with the ability to use modelMatrix...
    // - not sure how to make it aware of panelID (need to resolve how filtering works,
    //   currently I think it's slower than it should be)
    const chart = useChart();
    const cols = chart.config.param;
    // using both ref and state here so we can access the current value in the event handlers
    // can probably simplify this...
    const [start, setStartX] = useState<P>([0,0]);
    const [end, setEndX] = useState<P>([0,0]);
    const setStart = useCallback((p: P) => {
        startRef.current = p;
        setStartX(p);
    }, []);
    const setEnd = useCallback((p: P) => {
        endRef.current = p;
        setEndX(p);
    }, []);
    const startRef = useRef<P>([0,0]);
    const endRef = useRef<P>([0,0]);
    const uiElement = useRef<HTMLDivElement>(null);
    const scale = useRegionScale();//todo: this is a hack... should have a better way of reasoning about transforms...
    const updateRange = useCallback(() => {
        if (!rangeDimension) return;
        const s = startRef.current;
        const t = endRef.current;
        //need to convert from model coordinates to data coordinates...
        //seems like this may already be done somewhere?
        const range1 = [Math.min(s[0], t[0]), Math.max(s[0], t[0])]; //x range
        const range2 = [Math.min(s[1], t[1]), Math.max(s[1], t[1])]; //y range
        const args = { range1, range2 };
        rangeDimension.filter('filterSquare', cols, args);
        chart.resetButton.style.display = 'inline';
        (window as any).r = rangeDimension;
    }, [rangeDimension, cols, scale]);
    const unproject = useCallback((e: MouseEvent | React.MouseEvent) => {
        if (!uiElement.current) return [0,0] as P;
        const r = uiElement.current.getBoundingClientRect();
        const x = e.clientX - r.left;
        const y = e.clientY - r.top;
        const p = scatterplotLayer.unproject([x, y]) as P;
        //need to fix this now that we're using modelMatrix rather than scaling coordinates...
        return p.map(v => v*scale) as P;
    }, [scatterplotLayer, uiElement, scale]);
    const handleMouseMove = useCallback((e: MouseEvent) => {
        if (!toolActive) return;
        const p = unproject(e);
        setEnd(p);
    }, [toolActive]);
    const handleMouseUp = useCallback((e: MouseEvent) => {
        handleMouseMove(e);
        window.removeEventListener('mouseup', handleMouseUp);
        window.removeEventListener('mousemove', handleMouseMove);
        updateRange();
    }, [toolActive]);
    
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
        //consider using a ref to this element & updating the style directly...
        //to avoid needing to re-render the whole thing and keep both state and refs in sync...
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
    ref={uiElement}
    onMouseDown={(e) => {
        if (!toolActive) return;
        const p = unproject(e);
        setStart(p);
        setEnd(p);
        // setDragging(true); //dragging state is determined by whether the listeners are attached...
        window.addEventListener('mouseup', handleMouseUp);
        window.addEventListener('mousemove', handleMouseMove);
    }}
    />
    </>);
}


export default observer(function SelectionOverlay({scatterplotLayer} : {scatterplotLayer: ScatterplotLayer}) {
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
    const chart = useChart();
    const ds = useMemo(() => chart.dataStore, [chart]);
    const id = useChartID();
    const [rangeDimension, setRangeDimension] = useState<RangeDimension>(undefined);
    useEffect(() => {
        if (!ds) return;
        const rd = ds.getDimension('range_dimension');
        chart.removeFilter = () => {
            rd.removeFilter();
        }
        setRangeDimension(rd);

        return () => {
            chart.removeFilter = () => {};
            rd.destroy();
        }
    }, [ds]);


    // state: { selectedTool: 'rectangle' | 'circle' | 'polygon' | 'lasso' | 'magic wand' | 'none' }
    // interaction phases... maybe revert back to pan after making selection
    // - but there should be interaction with drag handles...
    // add or remove from selection...
    // -> transform into Deck coordinates...
    // later: selection layers...
    const drawRect = scatterplotLayer && scatterplotLayer.internalState; //<< I think this is causing unmounting issues...
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
                {drawRect && <RectangleEditor toolActive={selectedTool==='Rectangle'} scatterplotLayer={scatterplotLayer}
                rangeDimension={rangeDimension}
                />}
        </div>
        </>
    )
});