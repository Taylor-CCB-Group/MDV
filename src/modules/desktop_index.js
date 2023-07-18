import "./all_css"
import ChartManager from '../charts/ChartManager.js';
import { getArrayBufferDataLoader,getLocalCompressedBinaryDataLoader} from "../dataloaders/DataLoaders.js";

function _mdvInit(staticFolder){
    //get the configs for MDV
    getConfigs(staticFolder).then(resp=>{
        const config=resp.state;
        //is view in the URL
        const urlParams = new URLSearchParams(window.location.search);
        const view = urlParams.get("view");
        if (config.all_views && view && config.all_views.indexOf(view) !==-1){
            config["initial_view"]=view;
        }
        //data loaders depend on whether data is static or retrieved via API
        const dataLoader={
            function:staticFolder?getLocalCompressedBinaryDataLoader(resp.datasources,".")
                                    :getArrayBufferDataLoader("/get_data"),
            viewLoader:staticFolder?async (view)=> resp.views[view]
                                    :getView,
            rowDataLoader:staticFolder?loadRowDataStatic
                                    :loadRowData,
            binaryDataLoader:staticFolder?loadBinaryDataStatic
                                    :loadBinaryData
        };
        //lsiten to events
        const listener = (type,cm,data)=>{
            switch(type){
                case "state_saved":
                    getData("/save_state",data).then(resp=>{
                        if (resp.success){
                            cm.createInfoAlert("Data Saved",{duration:2000});
                            cm.setAllColumnsClean();
                        }
                        else{
                            cm.createInfoAlert("UnableToSaveData",{duration:3000,type:"danger"});
                        }
                    })	
                    break;
                case "view_loaded":
                    changeURLParam("view",cm.currentView);
                    break;
                }     
        };
        //create the app
        console.log(config);
        new ChartManager("holder",resp.datasources,dataLoader,config,listener);
    });
}
//only method required
window._mdvInit=_mdvInit;

//loads unstructured data for each row
async function loadRowData(datasource,index){
    return await getData("/get_row_data",{datasource,index})
}
async function loadRowDataStatic(datasource,index){
    const resp = await fetch(`./rowdata/${datasource}/${index}.json`);
    if (resp.status !=200){
        return null
    }
    return await resp.json()
}
//load view from API
async function getView(view){
    return await getData("/get_view",{view:view})
}

//load arbritray data
async function loadBinaryDataStatic(datasource,name){
    const resp = await fetch(`./binarydata/${datasource}/${name}.b`);
    return await resp.arrayBuffer();
}
async function loadBinaryData(datasource,name){
    return await getData("/get_binary_data",{datasource,name},"arraybuffer");
}

//get the configs wither from a folder or via a remote API
async function getConfigs(folder){
    let configs={};
    if (folder){
        let resp = await fetch(`./datasources.json`);
        configs.datasources = await resp.json();
        resp = await fetch(`./state.json`);
        configs.state = await resp.json();
        resp = await fetch(`./hyperion_config.json`);
        if (resp.status==200){
            configs.hyperion_config = await resp.json();
        }   
        resp = await fetch(`./views.json`);
        configs.views =  await resp.json();  
    }
    else{
        configs = await getData("/get_configs")
    }
    return configs;
}
//send json args and return json/array buffer response
async function getData(url,args,return_type="json"){
    const resp = await fetch(url,
    {
        method: "POST",
        body: JSON.stringify(args),
        headers: {
            "Accept": "application/json,text/plain,*/*",
            "Content-Type": "application/json"
        }
    });
    if (return_type==="json"){
        return await resp.json();
    }
    else{
        return await resp.arrayBuffer();
    }
    
}
//changes or adds a param to the browser address bar
function changeURLParam(param,value){
    const url = new URL(window.location);
    (url.searchParams.has(param) ? url.searchParams.set(param, value) : url.searchParams.append(param, value));
    url.search = url.searchParams;
    history.pushState({}, null, url);
}
