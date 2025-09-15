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
    return {
        data:{
            attributes:{
                getSV: svBuff,
                getColor: new Uint8Array(dataStore.size*3),
                getFilter: dataStore.filterArray,
            }
        },
        totalLength: cumSize
    };

}


const SVCircosPlotComponent = observer(()=>{
    const chart = useChart();
    const cb = chart.config.color_by;
    
    const {data,totalLength} = useMemo(()=>calculateSVs(chart.config.param,chart.dataStore),[]);
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
        const plotDiameter = (2 * 1000)+100;
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
    
        const layers = [
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