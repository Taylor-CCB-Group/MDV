import { ButtonGroup, IconButton, Tooltip } from "@mui/material"; //check tree-shaking...
import PanToolOutlinedIcon from '@mui/icons-material/PanToolOutlined';
import PhotoSizeSelectSmallOutlinedIcon from '@mui/icons-material/PhotoSizeSelectSmallOutlined';
import PolylineOutlinedIcon from '@mui/icons-material/PolylineOutlined';
import EditOutlinedIcon from '@mui/icons-material/EditOutlined';
import ControlCameraOutlinedIcon from '@mui/icons-material/ControlCameraOutlined';
import StraightenIcon from '@mui/icons-material/Straighten';
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useMetadata, useViewerStore } from "./avivatorish/state";
import { useFilteredIndices, useRegionScale, useScatterplotLayer } from "../scatter_state";
import { useChart } from "../context";
import { useMeasure, useRange } from "../spatial_context";
import RangeDimension from "../../datastore/RangeDimension";
import { observer } from "mobx-react-lite";
import type { VivMDVReact } from "./VivMDVReact";
import { runInAction } from "mobx";
import { useChartSize } from "../hooks";
import { sizeToMeters } from "./avivatorish/utils";

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
    'transform': {
        name: 'Transform',
        ToolIcon: ControlCameraOutlinedIcon
    },
    'measure': {
        name: 'Measure',
        ToolIcon: StraightenIcon
    },
} as const;

type Tool = typeof Tools[keyof typeof Tools]['name'];
const ToolArray = Object.values(Tools);
type P = [number, number];
type EditorProps = { toolActive?: boolean, rangeDimension: RangeDimension } & ReturnType<typeof useScatterplotLayer>;
function RectangleEditor({toolActive = false, scatterplotLayer, rangeDimension, unproject, currentLayerHasRendered} : EditorProps) {
    const chart = useChart() as VivMDVReact;
    const cols = chart.config.param;
    // using both ref and state here so we can access the current value in the event handlers
    // (without needing to recreate them every time the state changes)
    const { setStart, startRef, setEnd, endRef } = useRange();
    const { start, end } = useRange();
    // this is glitchy - not clear whether it was actually working correctly prior to local changes or not
    // ??? especially bad when measure tool is used first...
    const updateRange = useCallback(async () => {
        if (!rangeDimension) return;
        const s = startRef.current;
        const t = endRef.current;
        //need to convert from model coordinates to data coordinates...
        //seems like this may already be done somewhere?
        const range1 = [Math.min(s[0], t[0]), Math.max(s[0], t[0])]; //x range
        const range2 = [Math.min(s[1], t[1]), Math.max(s[1], t[1])]; //y range

        // this was not reflecting the background filter... so we select points that are not part of this region.
        // there needs to be more of a review of how our filters are evaluated...

        // the signature for this should be `(index: number) => boolean`
        // This 2D range predicate is responsible for the logic about which data to refer to.
        // columns argument to filter('filterPredicate', columns, args) is not used in the current implementation.
        // (unless there's some voodoo that signals to load column data)
        const data1 = rangeDimension.parent.columnIndex[cols[0]].data;
        const data2 = rangeDimension.parent.columnIndex[cols[1]].data;
        const predicate = (i: number) => {
            //filtered indices already include the filtering we do here...
            //we need a different strategy for this...
            // if (!indexSet.has(i)) return true;
            const v1 = data1[i];
            const v2 = data2[i];
            return !(v1 < range1[0] || v1 > range1[1] || v2 < range2[0] || v2 > range2[1] || isNaN(v1) || isNaN(v2))
        }
        const args = { range1, range2, predicate };

        //make zoom_on_filter only apply to other charts.
        const zoom_on_filter = chart.config.zoom_on_filter;
        runInAction(() => {
          chart.config.zoom_on_filter = false;
        });

        rangeDimension.filter('filterPredicate', cols, args);
        //although the filter is sync, the event that checks zoom_on_filter will happen later...
        //in particular, it runs in an effect that depends on async getFilteredIndices...
        //this is not a correct way to do this... also still needs testing with viewState link.
        //may be another place where react query would be useful?
        // await chart.dataStore.getFilteredIndices();
        setTimeout(() => {
            runInAction(() => {
              chart.config.zoom_on_filter = zoom_on_filter;
            });
        }, 500);
        chart.resetButton.style.display = 'inline';
        // (window as any).r = rangeDimension;
    }, [rangeDimension, cols, start, end]); //for some reason passing start & end here is improving the behaviour even though we use refs...

    const handleMouseMove = useCallback((e: MouseEvent) => {
        if (!toolActive) return;
        const p = unproject(e);
        setEnd(p);
    }, [toolActive, end]);
    const handleMouseUp = useCallback((e: MouseEvent) => {
        handleMouseMove(e);
        window.removeEventListener('mouseup', handleMouseUp);
        window.removeEventListener('mousemove', handleMouseMove);
        updateRange();
    }, [toolActive, updateRange]);

    if (!currentLayerHasRendered) return null; //if we pass this, I thought it meant we have internalState, but it seems not...
    if (!scatterplotLayer.internalState) return null;

    return (
    <>
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
        console.log('mouse down', p);
        setStart([p[0], p[1]]);
        setEnd([p[0], p[1]]);
        // setDragging(true); //dragging state is determined by whether the listeners are attached...
        window.addEventListener('mouseup', handleMouseUp);
        window.addEventListener('mousemove', handleMouseMove);
    }}
    />
    </>);
}

function MeasureTool({scatterplotLayer, unproject, toolActive} : EditorProps) {
    // click to set start, click to set end, draw line between them.
    // UX for tools is generally in need of work...
    const scale = useRegionScale();
    const metadata = useMetadata();
    const physicalSize = metadata?.Pixels?.PhysicalSizeX || 1;
    // const physicalUnits = metadata?.Pixels?.PhysicalSizeXUnit || 'unknown units';
    const [w, h] = useChartSize();
    const [hasStarted, setHasStarted] = useState(false);
    const { startPixels, setStart, endPixels, setEnd } = useMeasure();
    const handleMouseMove = useCallback((e: MouseEvent) => {
        setEnd(unproject(e));
    }, [unproject]);
    const handleMouseUp = useCallback((e: MouseEvent) => {
        handleMouseMove(e);
        window.removeEventListener('mouseup', handleMouseUp);
        window.removeEventListener('mousemove', handleMouseMove);
    }, [handleMouseMove]);
    if (!toolActive && !hasStarted) return null;
    // could we unproject into the image layer rather than scatterplotLayer?
    // const startPixels = unproject(start);
    // const endPixels = unproject(end);
    // todo make this a layer in deck.gl & hopefully avoid exceptions...
    const start = scatterplotLayer.project(startPixels);
    const end = scatterplotLayer.project(endPixels);
    const length = physicalSize * Math.sqrt((endPixels[0] - startPixels[0])**2 + (endPixels[1] - startPixels[1])**2) / scale;
    return (<>
    <div style={{
        position: 'absolute',
        width: '100%',
        height: '100px',
        display: hasStarted ? 'block' : 'none',
        top: end[1] + 10,
        left: end[0] + 10,
    }}>
        {(sizeToMeters(length, 'mm')).toFixed(2)}mm
    </div>
    <svg
    viewBox={`0 0 ${w} ${h}`}
    style={{
        position: 'absolute', top: 0, left: 0, width: '100%', height: '100%', backgroundColor: 'transparent',
        pointerEvents: toolActive ? 'auto' : 'none',
    }}
    onMouseDown={e => {
        if (!toolActive) return;
        const p = unproject([e.clientX, e.clientY]);
        setStart(p);
        setEnd(p);
        setHasStarted(true);
        window.addEventListener('mouseup', handleMouseUp);
        window.addEventListener('mousemove', handleMouseMove);        
    }}
    onScroll={e => {e.preventDefault();}}
    >
        {hasStarted && (<>
            <circle cx={start[0]} cy={start[1]} r={3} fill="white" />
            <circle cx={end[0]} cy={end[1]} r={3} fill="white" />
            <line x1={start[0]} y1={start[1]} x2={end[0]} y2={end[1]} stroke="white" strokeWidth={1.5}/>
        </>)}
    </svg>
    </>)
}

function TransformEditor({scatterplotLayer, modelMatrix, unproject} : EditorProps) {
    const pLastRef = useRef([NaN, NaN]);
    const handleMouseMove = useCallback((e: MouseEvent) => {
        const p = unproject(e);
        const pLast = pLastRef.current;
        const dx = p[0] - pLast[0];
        const dy = p[1] - pLast[1];
        modelMatrix.translate([dx, dy, 0]);
        //this may redraw nice & fast without incurring react overhead, but now we're left with other problems i.e. rectangle editor not updating...
        scatterplotLayer.setNeedsRedraw();
    }, []);
    const handleMouseUp = useCallback((e: MouseEvent) => {
        window.removeEventListener('mouseup', handleMouseUp);
        window.removeEventListener('mousemove', handleMouseMove);
    }, []);

    return (
        <>
        <div style={{
            position: 'absolute',
            width: '100%',
            height: '100%',
            top: 0,
            left: 0,
        }}
        onMouseDown={(e) => {
            const p = unproject(e);
            pLastRef.current = p;
            window.addEventListener('mouseup', handleMouseUp);
            window.addEventListener('mousemove', handleMouseMove);
        }}
        />
        </>
    )
}

export default observer(function SelectionOverlay(scatterProps : ReturnType<typeof useScatterplotLayer>) {
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
    const { rangeDimension } = useRange();
    // state: { selectedTool: 'rectangle' | 'circle' | 'polygon' | 'lasso' | 'magic wand' | 'none' }
    // interaction phases... maybe revert back to pan after making selection
    // - but there should be interaction with drag handles...
    // add or remove from selection...
    // -> transform into Deck coordinates...
    // later: selection layers...
    // glitches with settings changes are worse than the frame-behind overlay...
    // but if we render the overlay as a deck layer, we can avoid this... need to make a store for that.
    // if (!scatterProps.currentLayerHasRendered) return null;
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
            tabIndex={0}
            onMouseUp={(e) => {
                // setSelectedTool('Pan');
            }}
            onMouseDown={(e) => {
                if (selectedTool === 'Pan') {
                    // allow the event to propagate to the deck canvas below
                    return;
                }
                e.currentTarget.focus();
            }}
            onBlurCapture={(e) => {console.log('blur capture');}}
            onBlur={(e) => {console.log('blur');}}
            onKeyDownCapture={(e) => {
                if (e.key === 'Escape' || e.key === ' ') {
                    setSelectedTool('Pan');
                }
                if (e.key === 'r') {
                    setSelectedTool('Rectangle');
                }
                if (e.key === 't') {
                    setSelectedTool('Transform');
                }
                if (e.key === 'm') {
                    setSelectedTool('Measure');
                }
            }}
            >
                <RectangleEditor toolActive={selectedTool==='Rectangle'} {...scatterProps}
                rangeDimension={rangeDimension} />
                {selectedTool === 'Transform'  && <TransformEditor {...scatterProps} rangeDimension={rangeDimension}/>}
                {scatterProps.currentLayerHasRendered && <MeasureTool {...scatterProps} rangeDimension={rangeDimension} toolActive={selectedTool === 'Measure'} />}
        </div>
        </>
    )
});
