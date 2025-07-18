

import { getChart } from "./chartTypes";
import {useMemo,useEffect,useState,useRef} from "react";
import { v4 as uuidv4 } from 'uuid';
import { Instance } from "mobx-state-tree";
import { link } from "node:fs";
import { columnMatchesType } from "@/lib/columnTypeHelpers";
type Filter = {
  timeStamp: number;
  dimension: any; // Replace 'any' with a more specific type if possible
};

type MDVChartProps = {
  chartInstance: Instance<any>;
  datastore: any;
  filter: Filter
};

//temporary function to get the datastore from the chartManager
function getDataStore(datastore: string) {
    const cm = window.mdv.chartManager;
    return cm.getDataSource(datastore);
}
function getColumn(datastore: string, field: string) {
    return new Promise((resolve, reject) => {
        const cm = window.mdv.chartManager;
        cm._getColumnsThen(datastore, [field],() => {
            resolve(cm.getDataSource(datastore).columnIndex[field]);
        })
                        
    });
}


const useColumn= (column: Instance<any>) => {
    const [data, setData] = useState<any>(null);
    const linkRef = useRef<{ id: string; ds: any } | null>(null);
    

    useEffect(() => {
        getColumn(column.datastore,column.field).then((result: any) => {
            setData(result);
        });
    }, [column.field,column.datastore]);
    useEffect(()=>{
        
        if (linkRef.current) {  
            linkRef.current.ds.removeListener(linkRef.current.id);
        }
        if (column.link){
            const id =  uuidv4();
            const dst  = getDataStore(column.link.datastore);
            const ds =    getDataStore(column.datastore);
            dst.addListener(id, async (type: string, data: any) => {
                if (type === "data_highlighted") {
                    console.log(data);
                }
                const n_col = ds.links?.[column.link.datastore]?.rows_as_columns?.name_column;
                if (typeof n_col === "string") {
                    await getColumn(column.link.datastore, n_col);
                    const i = data.indexes[0];
                    const n  = dst.getRowText(i, n_col);
                    const g= column.link.subgroup;
                    column.setField(`${g}|${n} (${g})|${i}`);

                }
                console.log("here");
            })
            linkRef.current = {
                id,
                ds: dst
            };
        }
        return ()=>{
            if (linkRef.current) {
                linkRef.current.ds.removeListener(linkRef.current.id);
            }
        }
    },[column.link])
    return data;
}

const useFilter = (datastore: any) => {

    const [filter, setFilter] = useState<Filter | {}>({});
    useEffect(() => {
        const id  = uuidv4()
        datastore.addListener(id,(type: string,data: any,rows:number)=>{
            if (type==="filtered"){
                //if a dimension is doing the filtering then this will
                //be in the data else undefined
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



export default function MDVChart({ config, datastore }: { config: any, datastore: any }): React.ReactElement {
    // Retrieve the chart type entry using the config type
    // and create an instance of the chart model with the provided config.
    const filter = useFilter(datastore);
    const { chartInstance, ChartComponent } = useMemo(() => {
        const chartEntry = getChart(config.type);
        if (!chartEntry) {
            return { chartInstance: null, ChartComponent: null };
        }
        const { chartModel, ChartComponent } = chartEntry;
        return { chartInstance: chartModel.create(config), ChartComponent };
    }, []);
  return (
    // menu bar /title bar goes here
    <>
        {ChartComponent
        ? <ChartComponent chartInstance={chartInstance} datastore={datastore} filter={filter}/>
        : <h2>Chart type not found: {config.type}</h2>}
    </>
  );
}

export { useFilter, useColumn};
export type { MDVChartProps, Filter };