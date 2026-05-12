import {observer} from "mobx-react-lite";
import SVLayer from "./SVLayer.js";
import {useState} from "react";
import DeckGL from 'deck.gl';
import { OrthographicView } from '@deck.gl/core';
import { v4 as uuidv4 } from 'uuid';
import {useEffect,useMemo} from "react";
import {useOuterContainer} from "../../react/screen_state";

import {useChart} from "../../react/context";
import { useChartSize } from "@/react/hooks.js";
import { PathLayer, TextLayer } from "@deck.gl/layers";

const CHROMOSOME_COLORS = [
    [31, 119, 180],
    [255, 127, 14],
    [44, 160, 44],
    [214, 39, 40],
    [148, 103, 189],
    [140, 86, 75],
    [227, 119, 194],
    [127, 127, 127],
    [188, 189, 34],
    [23, 190, 207],
];


//don't want to use filtredindxes as the filter buffer is used directly]
//in wbgl , just need to know when filtering happans
const useFilter = (datastore) => {
    const [filter, setFilter] = useState({});
    useEffect(() => {
        const id  = uuidv4()
        datastore.addListener(id,(type, data)=>{
            if (type==="filtered"){
                //if a dimension is doing the filtering then this will
                //be in the data object  else undefined
                setFilter({
                    timeStamp: Date.now(),
                    dimension:data
                })
            }
        })
        return () =>{
            datastore.removeListener(id);
        }
    },[])
  return filter;
}

//convert svs in into a flat array of start, length, type
//for webgl - should move this to a web worker if it gets too slow
function calculateSVs(params,dataStore){
    // get the chromosome in proper order 
    const chrs = Object.keys(dataStore.genome.chromosomes).sort((a, b) => 
        a.localeCompare(b, undefined, { numeric: true, sensitivity: 'base' })
    );
    const chrCumLen = {};
    let cumSize = 0;
    for (let i = 0; i < chrs.length; i++) {
        const chr = chrs[i];
        chrCumLen[chr] = cumSize;
        cumSize += dataStore.genome.chromosomes[chr];
    }
    const sv_map={
            "TRA":0,
            "BND":0,
            "INV":1,
            "DEL":2,
            "INS":3,
        " DUP":4,
    }

    const p =params
    const pos1  =  dataStore.getRawColumn(p[0]);
    const chr1  =  dataStore.getRawColumn(p[1]);
    const pos2  =  dataStore.getRawColumn(p[2]);
    const chr2  =  dataStore.getRawColumn(p[3]);
    const svtype=  dataStore.getRawColumn(p[4]);
    const chr1Vals= dataStore.getColumnValues(p[1]);
    const chr2Vals= dataStore.getColumnValues(p[3]);
    const svtypeVals= dataStore.getColumnValues(p[4]);
    const svBuff = new Uint32Array(dataStore.size*3);
    for (let i = 0; i < dataStore.size; i++) {
        const idx= i*3;
        const sv_start = pos1[i] + chrCumLen[chr1Vals[chr1[i]]];
        const sv_end   = pos2[i] + chrCumLen[chr2Vals[chr2[i]]];
        let t = sv_map[svtypeVals[svtype[i]]];
        if (t===undefined){
            t=5;
        }
        svBuff[idx] = sv_start;
        svBuff[idx+1] = sv_end -sv_start;
        svBuff[idx+2] = t;
    }
    const chromosomeSegments = chrs.map((chr, index) => ({
        chr,
        start: chrCumLen[chr],
        length: dataStore.genome.chromosomes[chr],
        index
    }));

    return {
        data:{
            attributes:{
                getSV: svBuff,
                getColor: new Uint8Array(dataStore.size*3),
                getFilter: dataStore.filterArray,
            }
        },
        totalLength: cumSize,
        chromosomeSegments
    };

}

function getChromosomePath(start, length, totalLength, radius, steps = 48) {
    const startAngle = (start * 2 * Math.PI) / totalLength - Math.PI / 2;
    const endAngle = ((start + length) * 2 * Math.PI) / totalLength - Math.PI / 2;
    const stepCount = Math.max(2, steps);
    const points = [];
    for (let i = 0; i <= stepCount; i++) {
        const t = i / stepCount;
        const angle = startAngle + (endAngle - startAngle) * t;
        points.push([Math.cos(angle) * radius, Math.sin(angle) * radius]);
    }
    return points;
}

function pointAtGenomePosition(position, totalLength, radius) {
    const angle = (position * 2 * Math.PI) / totalLength - Math.PI / 2;
    return [Math.cos(angle) * radius, Math.sin(angle) * radius];
}

function getChromosomeLabelData(chromosomeSegments, totalLength, radius) {
    return chromosomeSegments.map((segment) => {
        const mid = segment.start + segment.length / 2;
        return {
            text: segment.chr,
            position: pointAtGenomePosition(mid, totalLength, radius),
        };
    });
}

function formatTickLabel(bp) {
    if (bp >= 1000000) return `${Math.round(bp / 1000000)}Mb`;
    if (bp >= 1000) return `${Math.round(bp / 1000)}kb`;
    return `${bp}bp`;
}

function getTickData(chromosomeSegments, totalLength, tickStep = 25000000, majorEvery = 2) {
    const outerTickLines = [];
    const outerTickLabels = [];

    chromosomeSegments.forEach((segment) => {
        const ticks = Math.floor(segment.length / tickStep);
        for (let i = 0; i <= ticks; i++) {
            const offset = i * tickStep;
            if (offset > segment.length) break;
            const genomePos = segment.start + offset;
            const pOuter1 = pointAtGenomePosition(genomePos, totalLength, 1052);
            const isMajor = i % majorEvery === 0;
            const pOuter2 = pointAtGenomePosition(genomePos, totalLength, isMajor ? 1080 : 1066);
            outerTickLines.push({
                path: [pOuter1, pOuter2],
                major: isMajor,
            });

            if (isMajor && i > 0) {
                const outerLabelPos = pointAtGenomePosition(genomePos, totalLength, 1095);
                outerTickLabels.push({
                    text: formatTickLabel(offset),
                    position: outerLabelPos,
                });
            }
        }
    });
    return { outerTickLines, outerTickLabels };
}

function getDisplaySpacing(zoomLevel) {
    if (zoomLevel < -0.4) return { tickStep: 100000000, majorEvery: 1 };
    if (zoomLevel < 0.3) return { tickStep: 50000000, majorEvery: 1 };
    if (zoomLevel < 1.0) return { tickStep: 25000000, majorEvery: 2 };
    return { tickStep: 10000000, majorEvery: 2 };
}


const SVCircosPlotComponent = observer(()=>{
    const chart = useChart();
    const cb = chart.config.color_by;
    
    const {data,totalLength,chromosomeSegments} = useMemo(()=>calculateSVs(chart.config.param,chart.dataStore),[]);
    const doc = useOuterContainer();
    const [width,height] = useChartSize();
    

    //seems to be called twice when cb changes
    //also will not update when other aspects of color function change
    useMemo(()=>{
        if (cb){
            const colorFunction = chart.getColorFunction(cb,true);
            const buf = data.attributes.getColor;
            for (let i = 0; i < chart.dataStore.size; i++) {
                const color = colorFunction(i);
                const idx = i*3;
                buf[idx] = color[0];
                buf[idx+1] = color[1];
                buf[idx+2] = color[2];
            }
        }
    },[cb])
    const [viewState, setViewState] = useState({
          target: [0, 0],
          zoom: [0, 0] // [xZoom, yZoom]
    });

    useEffect(() => {
        // Fit the plot to the container
        const plotDiameter = (2 * 1120)+100;
        const minDim = Math.min(width, height);
        const zoom = Math.log2(minDim / plotDiameter);

        setViewState({
            target: [0, 0], // center of plot
            zoom: [zoom, zoom], // fit both axes
        });
    }, [width, height]);

    const [hoveredIndex, setHoveredIndex] = useState(-1);
    const filter = useFilter(chart.dataStore);

    const handleHover = (info) => {
        setHoveredIndex(info && info.index >= 0 ? info.index : -1);
    };

    const clicked = (info, event) => {
      if (hoveredIndex !==-1){
        chart.dataStore.dataHighlighted([hoveredIndex],chart);
      } 
    }
    
        const zoomLevel = Array.isArray(viewState.zoom) ? viewState.zoom[0] : 0;
        const { tickStep, majorEvery } = getDisplaySpacing(zoomLevel);
        const showTicks = zoomLevel > -0.9;
        const showTickLabels = showTicks;

        const chromosomeLayer = new PathLayer({
            id: "chromosome_layer",
            data: chromosomeSegments,
            getPath: (d) => getChromosomePath(d.start, d.length, totalLength, 1045),
            getColor: (d) => CHROMOSOME_COLORS[d.index % CHROMOSOME_COLORS.length],
            getWidth: 8,
            widthUnits: "pixels",
            widthMinPixels: 2,
            pickable: false,
            capRounded: false,
            jointRounded: false
        });

        const chromosomeLabels = getChromosomeLabelData(chromosomeSegments, totalLength, 1120);
        const {
            outerTickLines,
            outerTickLabels
        } = getTickData(chromosomeSegments, totalLength, tickStep, majorEvery);

        const chromosomeLabelLayer = new TextLayer({
            id: "chromosome_label_layer",
            data: chromosomeLabels,
            getPosition: (d) => d.position,
            getText: (d) => d.text,
            getColor: [30, 30, 30, 255],
            getSize: 14,
            sizeUnits: "pixels",
            sizeMinPixels: 10,
            getTextAnchor: "middle",
            getAlignmentBaseline: "center",
            pickable: false
        });

        const tickLayer = new PathLayer({
            id: "chromosome_tick_layer",
            data: outerTickLines,
            getPath: (d) => d.path,
            getColor: (d) => (d.major ? [70, 70, 70, 255] : [120, 120, 120, 220]),
            getWidth: (d) => (d.major ? 2.0 : 1.0),
            widthUnits: "pixels",
            widthMinPixels: 1,
            visible: showTicks,
            pickable: false
        });

        const outerTickLabelLayer = new TextLayer({
            id: "chromosome_tick_outer_label_layer",
            data: outerTickLabels,
            getPosition: (d) => d.position,
            getText: (d) => d.text,
            getColor: [80, 80, 80, 255],
            getSize: 10,
            sizeUnits: "pixels",
            sizeMinPixels: 8,
            getTextAnchor: "middle",
            getAlignmentBaseline: "center",
            visible: showTickLabels,
            pickable: false
        });

        const layers = [
        chromosomeLayer,
        tickLayer,
        outerTickLabelLayer,
        chromosomeLabelLayer,
        new SVLayer({
            id: 'sv_layer',
            data,
            pickable:true,
            highlightColor: [0,0, 0, 255],
            highlightedObjectIndex: hoveredIndex,
            radius: 1000,
            totalLength,
            container:doc.body,
            center: [0, 0],
            rotation:0,
            updateTriggers:{
              getFilter:[filter.timeStamp],
              getColor:[cb]
            }    
        })
        ];

    const handleViewStateChange = ({viewState: nextViewState}) => {
        setViewState(prev => ({
          ...nextViewState,
        }));
      };
    return(
        <DeckGL
            layers={layers}
            viewState={viewState}
            onViewStateChange={handleViewStateChange}
            onClick={clicked}
            onHover={handleHover}
            getCursor={() => hoveredIndex >= 0 ? 'pointer' : 'default'}
            container={document.body}
            views={new OrthographicView()}
            controller={true}
        />  
    )
})



export default SVCircosPlotComponent;
