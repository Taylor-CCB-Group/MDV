import "../modules/all_css.js"
import ChartManager from '../charts/ChartManager.js';
import {getLocalCompressedBinaryDataLoader} from "../dataloaders/DataLoaders.js";

async function init(){
    let r =  await fetch("datasources.json");
    const ds = await r.json();
    r= await fetch("views.json");
    const vs= await r.json();
    r= await fetch("state.json");
    const cf = await r.json();
    const dataLoader={
        function:getLocalCompressedBinaryDataLoader(ds,"."),
        viewLoader:async v=> vs[v],
        rowDataLoader:async (datasource,index)=>{
            const resp = await fetch(`./rowdata/${datasource}/${index}.json`);
            if (resp.status !=200){
                return null;
            }
            return await resp.json()
        } ,
        binaryDataLoader:async (datasource,name)=>{
            const resp = await fetch(`./binarydata/${datasource}/${name}.b`);
            return resp.arrayBuffer();
        }

    }
    new ChartManager("holder",ds,dataLoader,cf);
}

init();