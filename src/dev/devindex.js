import 'nouislider/dist/nouislider.min.css'
import "microtip/microtip.css";
import "../css/fontawesome-5.15.3/all.css";
import "../utilities/css/ContextMenu.css";
import "../charts/css/charts.css";
import "../webgl/css/wgl2di.css";
import "../table/css/slickgrid.css";
import "../browser/css/browser.css";
import "../browser/bam_track.js";
import "../browser/BamCoverageTrack.js";
import "../charts/GenomeBrowser.js";
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
        } 
    }
    new ChartManager("holder",ds,dataLoader,cf);
}


init();
